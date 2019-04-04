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

package Bio::EnsEMBL::Variation::Utils::RunDbNSFPUtils;

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

  my ($working_dir, $dbnsfp_file, $dbnsfp_version, $assembly) = rearrange([qw(WORKING_DIR DBNSFP_FILE DBNSFP_VERSION ASSEMBLY)], @_);
  $self->{'working_dir'} = $working_dir;
  $self->{'dbnsfp_file'} = $dbnsfp_file;
  $self->{'dbnsfp_version'} = $dbnsfp_version;
  $self->{'assembly'} = $assembly;
  if (! grep {$_ eq $assembly} ('GRCh37', 'GRCh38')) {
    die "Assembly $assembly is not supported.";
  }
  if (! grep {$_ eq $dbnsfp_version} ('3.5a')) {
    die "dbNSFP version $dbnsfp_version is not supported.";
  }

  return $self;
}

my $predictions = {
  dbnsfp_meta_lr => {
    T => 'tolerated',
    D => 'damaging',
  },
  dbnsfp_mutation_assessor => {
    H => 'high',
    M => 'medium',
    L => 'low',
    N => 'neutral',
  }
};

my $column_names = {
  '3.5a' => {
    assembly_unspecific => {
      chr => '#chr',
      ref => 'ref',
      refcodon => 'refcodon',
      alt => 'alt',
      aaalt => 'aaalt',
      aaref => 'aaref',
      revel_score => 'REVEL_score',
      meta_lr_score => 'MetaLR_score',
      meta_lr_pred => 'MetaLR_pred',
      mutation_assessor_score => 'MutationAssessor_score_rankscore',
      mutation_assessor_pred => 'MutationAssessor_pred',
    },
    'assembly_specific' => {
      'GRCh37' => {
        pos => 'pos(1-based)'
      },
      'GRCh38' => {
        pos => 'pos(1-based)'
      },
    },
  },
  '4.0b1' => {

  }
};

sub run {
  my $self = shift;
  my $translation_md5 = shift;

  my $translation_stable_id = $self->get_stable_id_for_md5($translation_md5);
  my $translation = $self->get_translation($translation_stable_id);
  my $translation_seq = $translation->seq;
  my $transcript = $translation->transcript;
  my $reverse = $transcript->strand < 0;
  my $transcript_stable_id = $transcript->stable_id;

  $self->init_protein_matrix($translation, $translation_md5);

  my @amino_acids = ();
  my @all_triplets = @{$self->get_triplets($translation_stable_id)};

  $self->init_header;
  my $obj = Bio::DB::HTS::Tabix->new(filename => $self->dbnsfp_file);

  my $codonTable = Bio::Tools::CodonTable->new();

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
      my $iter = $obj->query("$chrom:$triplet_start-$triplet_end");
      while (my $line = $iter->next) {
        my $data = $self->get_dbNSFP_row($line);
        my $chr = $data->{'chr'};
        my $pos = $data->{'pos'};
        my $ref = $data->{'ref'};
        my $refcodon = $data->{'refcodon'};
        my $alt = $data->{'alt'};
        my $aaalt = $data->{'aaalt'};
        my $aaref = $data->{'aaref'};
        next if ($alt eq $ref);

        my $nucleotide_position = ($reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $mutated_triplet =  $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $codonTable->translate($mutated_triplet);

        next if ($aaalt ne $mutated_aa);
        if ($data->{revel_score} ne '.') {
          my $prediction = ($data->{revel_score} >= $REVEL_CUTOFF) ? 'likely disease causing' : 'likely benign';
          $self->add_prediction($i, $mutated_aa, 'dbnsfp_revel', $data->{revel_score}, get_revel_prediction($prediction));
        }
        if ($data->{meta_lr_score} ne '.') {
          my $prediction = $predictions->{dbnsfp_meta_lr}->{$data->{meta_lr_pred}};
          $self->add_prediction($i, $mutated_aa, 'dbnsfp_meta_lr', $data->{meta_lr_score}, $prediction);
        }
        if ($data->{mutation_assessor_score} ne '.') {
          my $prediction = $predictions->{dbnsfp_mutation_assessor}->{$data->{mutation_assessor_pred}};
          $self->add_prediction($i, $mutated_aa, 'dbnsfp_mutation_assessor', $data->{mutation_assessor_score}, $prediction);
        }
      }    
    } # end foreach coord
  } # end foreach triplet
  if ($self->{'pipeline_mode'}) {
    if ($translation_seq ne join('', @amino_acids)) {
      my $fh = FileHandle->new($self->working_dir. "/$translation_stable_id", 'w');
      print $fh "$transcript_stable_id\n$translation_seq\n";
      print $fh join('', @amino_acids), "\n";
      $fh->close;
    }
  }

  $self->store_protein_matrix($translation_stable_id, $translation_md5);
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

sub get_dbNSFP_row {
  my $self = shift;
  my $line = shift;
  $line =~ s/\r$//g;
  my @split = split /\t/, $line;
  my $header = $self->header;
  my $assembly = $self->assembly;
  my $dbnsfp_version = $self->dbnsfp_version;
  my %raw_data = map {$header->[$_] => $split[$_]} (0..(scalar @{$header} - 1));
  my $data = {};
  my $assembly_unspecific = $column_names->{$dbnsfp_version}->{assembly_unspecific};
  foreach my $column_name (keys %{$assembly_unspecific}) {
    $data->{$column_name} = $raw_data{$assembly_unspecific->{$column_name}};
  }
  my $assembly_specific =  $column_names->{$dbnsfp_version}->{assembly_specific}->{$assembly}; 
  foreach my $column_name (keys %{$assembly_specific}) {
    $data->{$column_name} = $raw_data{$assembly_specific->{$column_name}};
  }
  return $data;
}

sub working_dir {
  my $self = shift;
  return $self->{'working_dir'};
}

sub dbnsfp_file {
  my $self = shift;
  return $self->{'dbnsfp_file'};
}

sub dbnsfp_version {
  my $self = shift;
  return $self->{'dbnsfp_version'};
}

sub assembly {
  my $self = shift;
  return $self->{'assembly'};
}

sub header {
  my $self = shift;
  return $self->{'header'} = shift if(@_);
  return $self->{'header'};
}

sub init_header {
  my $self = shift;
  my $header;
  open HEAD, "tabix -fh $self->dbnsfp_file 1:1-1 2>&1 | ";
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
  foreach my $analysis (qw/dbnsfp_revel dbnsfp_meta_lr dbnsfp_mutation_assessor/) {
    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
      -analysis       => $analysis,
      -peptide_length   => $translation->length,
      -translation_md5  => $translation_md5,
    );
    $pred_matrices->{$analysis} = $pred_matrix;
  }

  $self->{pred_matrices} = $pred_matrices;
}

sub store_protein_matrix {
  my $self = shift;
  my $translation_stable_id = shift;
  my $translation_md5 = shift;
  my $pred_matrices = $self->{pred_matrices};

  my $vdba = $self->get_species_adaptor('variation');
  my $pfpma = $vdba->get_ProteinFunctionPredictionMatrixAdaptor or die "Failed to get matrix adaptor";

  foreach my $analysis (keys %$pred_matrices) {
    my $pred_matrix = $pred_matrices->{$analysis};
    if ($self->{results_available}->{$analysis}) {
      $pfpma->store($pred_matrix);

      if ($self->{'debug_mode'}) {
        my $fh = FileHandle->new($self->working_dir. "/$analysis\_$translation_stable_id", 'w');
        my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, $translation_md5); 
        my $debug_data = $self->{debug_data};
        foreach my $i (keys %{$debug_data->{$analysis}}) {
          foreach my $aa (keys %{$debug_data->{$analysis}->{$i}}) {
            next if ($aa eq '*');
            foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
              my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
              print $fh join(' ', $analysis, $i, $aa, $prediction, $debug_data->{$analysis}->{$i}->{$aa}->{$prediction}, $new_pred, $new_score), "\n";
            }
          }
        }
        $fh->close;
      }
    }
  }
} 

1;
