=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut



=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 DESCRIPTION

This module contains functions used in the variant quality control process. 

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::RunCADDUtils;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionUtils;
use Bio::DB::HTS::Tabix;
use Bio::Tools::CodonTable;
use File::Path qw(make_path remove_tree);
use FileHandle;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);
our @ISA = ('Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionUtils');

my $REVEL_CUTOFF = 0.5;
my $LOW_QUALITY = 0;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($working_dir, $cadd_file, $dbnsfp_version, $assembly, $pipeline_mode, $debug_mode) = rearrange([qw(WORKING_DIR CADD_FILE CADD_VERSION ASSEMBLY PIPELINE_MODE DEBUG_MODE)], @_);
  $self->{'working_dir'} = $working_dir;
  $self->{'cadd_file'} = $dbnsfp_file;
  $self->{'cadd_version'} = $dbnsfp_version;
  $self->{'assembly'} = $assembly;
  if (! grep {$_ eq $assembly} ('GRCh37', 'GRCh38')) {
    die "Assembly $assembly is not supported.";
  }

  $self->{'pipeline_mode'} = (!defined $pipeline_mode) ? 1 : $pipeline_mode; # default set to 1
  $self->{'debug_mode'} = $debug_mode;

  return $self;
}

sub run {
  my $self = shift;
  my $translation_md5 = shift;

  my $translation_stable_id = $self->get_stable_id_for_md5($translation_md5);
  my $translation = $self->get_translation($translation_stable_id);
  my $translation_seq = $translation->seq;
  my $transcript = $translation->transcript;
  $self->reverse($transcript->strand < 0);
  my $transcript_stable_id = $transcript->stable_id;

  $self->init_protein_matrix($translation, $translation_md5);

  $self->init_header;

  my $all_triplets = $self->get_triplets($translation_stable_id);

  $self->load_predictions_for_triplets($all_triplets);

  if ($self->{'pipeline_mode'}) {
    if ($translation_seq ne join('', @{$self->amino_acids})) {
      my $fh = FileHandle->new($self->working_dir. "/$translation_stable_id", 'w');
      print $fh "$transcript_stable_id\n$translation_seq\n";
      print $fh join('', @{$self->amino_acids}), "\n";
      $fh->close;
    }
    $self->store_protein_matrix($translation_stable_id, $translation_md5);
  }
}

sub reverse {
  my $self = shift;
  return $self->{'reverse'} = shift if(@_);
  return $self->{'reverse'};
}

sub amino_acids {
  my $self = shift;
  my $aa = shift;
  if (defined $aa) {
    push @{$self->{'amino_acids'}}, $aa; 
  } else {
    return $self->{'amino_acids'};
  }
}

sub load_predictions_for_triplets {
  my $self = shift;
  my $triplets = shift; 
  foreach my $entry (@all_triplets) {
    my $aa = $entry->{aa};
    push @amino_acids, $aa;
    next if $aa eq 'X';
    my @coords = @{$entry->{coords}};
    my $chrom = $entry->{chrom};
    my $triplet_seq = $entry->{triplet_seq};
    my $i = $entry->{aa_position};
    my $new_triplets = $entry->{new_triplets};
    foreach my $coord (@coords) {
      my $triplet_start = $coord->[0];
      my $triplet_end = $coord->[1];
      my $iter = $self->get_tabix_iterator($chrom, $triplet_start, $triplet_end);
      while (my $line = $iter->next) {
        my $data = $self->get_CADD_row($line);
        my $chr = $data->{'#Chr'};
        my $pos = $data->{'Pos'};
        my $ref = $data->{'Ref'};
        my $alt = $data->{'Alt'};
        my $cadd_phred = $data{'PHRED'};
        next if ($alt eq $ref);
        my $nucleotide_position = ($reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $mutated_triplet = $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $codonTable->translate($mutated_triplet);
        $self->add_predictions($data, $i, $mutated_aa);
      }
    }
  }
} 

sub parser {
  my $self = shift;
  if (!defined $self->{'parser'}) {
    my $cadd_file = $self->cadd_file;
    $self->{'parser'} = Bio::DB::HTS::Tabix->new(filename => $cadd_file);
  }
  return $self->{'parser'};
}

sub codon_table {
  my $self = shift;
  if (!defined $self->{'codon_table'}) {
    $self->{'codon_table'} = Bio::Tools::CodonTable->new();
  }
  return $self->{'codon_table'};
}

sub get_tabix_iterator {
  my ($self, $chrom, $triplet_start, $triplet_end) = @_;
  return $self->parser->query("$chrom:$triplet_start-$triplet_end");
}

sub add_predictions {
  my ($self, $data, $i, $mutated_aa) = @_;
  if ($data->{'PHRED'} ne '.') {
    my $prediction = ($data->{'PHRED'} >= $CADD_CUTOFF) ? 'likely deleterious' : 'likely benign';
    $self->add_prediction($i, $mutated_add, 'cadd', $data->{'PHRED'}, $prediction);
  }
}

sub add_prediction {
  my ($self, $i, $mutated_aa, $predictor, $score, $prediction) = @_; 
  $self->{pred_matrices}->{$predictor}->add_prediction(
    $i,
    $mutated_aa,
    $prediction,
    $score,
    $LOW_QUALITY,
  );
 
  $self->{results_available}->{$predictor} = 1;
  $self->{debug_data}->{$predictor}->{$i}->{$mutated_aa}->{$prediction} = $score if ($self->{'debug_mode'});
}

sub get_CADD_row {
  my $self = shift;
  my $line = shift;
  $line =~ s/\r$//g;
  my @split = split /\t/, $line;
  my $header = $self->header;
  my %data = map {$header->[$_] => $split[$_]} (0..(scalar @{$header} - 1));
  return \%data;
}

sub working_dir {
  my $self = shift;
  return $self->{'working_dir'};
}

sub cadd_file {
  my $self = shift;
  return $self->{'cadd_file'};
}

sub cadd_version {
  my $self = shift;
  return $self->{'cadd_version'};
}

sub assembly {
  my $self = shift;
  return $self->{'assembly'};
}

sub header {
  my $self = shift;
  return $self->{'header'} = shift if(@_);
  $self->init_header if (!defined $self->{'header'});
  return $self->{'header'};
}

sub init_header {
  my $self = shift;
  my $header;
  my $cadd_file = $self->cadd_file;
  open HEAD, "tabix -fh $cadd_file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $header = [split];
  }
  close HEAD;

  $self->header($header);
}

sub init_protein_matrix {
  my $self = shift;
  my $translation = shift;
  my $translation_md5 = shift;
  my $pred_matrices = {};
  foreach my $analysis (qw/cadd/) {
    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
      -analysis       => $analysis,
      -peptide_length   => $translation->length,
      -translation_md5  => $translation_md5,
    );
    $pred_matrices->{$analysis} = $pred_matrix;
  }

  $self->{pred_matrices} = $pred_matrices;
}

1;
