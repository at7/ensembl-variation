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

package Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotation;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::Tools::CodonTable;
use Bio::DB::HTS::Tabix;
use Bio::EnsEMBL::Variation::Utils::BaseDatabaseUtils;
use FileHandle;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);

our @ISA = ('Bio::EnsEMBL::Variation::Utils::BaseDatabaseUtils');

my $LOW_QUALITY = 0;

=head2 new

  Arg [-working_dir] :
    string - location of the working directory when running the pipeline. The directory is used for storing debug information.
  Arg [-annotation_file] :
    string - location of dbNSFP or CADD file
  Arg [-annotation_file_version] :
    string - version of annotation file
  Arg [-assembly] :
    string - assembly version e.g. GRCh37 or GRCh38
  Arg [-pipeline_mode] :
    boolean - If set to 1 run in pipeline mode and store results in the database. The default value is 1.
  Arg [-debug_mode] :
    boolean - If set to 1 write debug information to the working directory. The default value is 0.
  Arg [-write_mode] :
    boolean - If set to 1 write error file which reports missmatch between protein sequence and translated protein sequence as a result of annotation.

  Example    : 
  my $cadd = Bio::EnsEMBL::Variation::Utils::RunCADDAnnotationUtils->new(
    -species => 'Homo_sapiens',
    -annotation_file => $dir . '/testdata/cadd_v1.3_grch37.txt.gz',
    -assembly => 'GRCh37',
    -annotation_file_version => 'v1.3',
    -pipeline_mode => 0,
    -debug_mode => 1,
  );

  Description: Constructor. Instantiates a new BaseProteinFunctionAnnotation object.
  Returntype : BaseProteinFunctionAnnotation 
  Exceptions : none
  Caller     : general
  Status     : Stable
=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  my ($working_dir, $annotation_file, $annotation_file_version, $assembly, $pipeline_mode, $debug_mode, $write_mode) = rearrange([qw(WORKING_DIR ANNOTATION_FILE ANNOTATION_FILE_VERSION ASSEMBLY PIPELINE_MODE DEBUG_MODE WRITE_MODE)], @_);
  $self->{'working_dir'} = $working_dir;
  $self->{'annotation_file'} = $annotation_file;
  $self->{'annotation_file_version'} = $annotation_file_version;
  $self->{'assembly'} = $assembly;
  $self->{'pipeline_mode'} = (!defined $pipeline_mode) ? 1 : $pipeline_mode; # default set to 1
  $self->{'write_mode'} = (!defined $write_mode) ? 1 : $write_mode; # default set to 1
  $self->{'debug_mode'} = $debug_mode;

  if (! grep {$_ eq $self->assembly} ('GRCh37', 'GRCh38')) {
    die "Assembly $assembly is not supported.";
  }

  return $self;
}

=head2 run
  Arg [1]    : String $translation_md5  
  Example    : $->run('');
  Description: Runs protein function prediction annotations for a protein translation.
  Returntype : none
  Exceptions : throws on missing argument
               throws on undefined $translation_stable_id
  Caller     : general
  Status     : At Risk
=cut

sub run {
  my $self = shift;
  my $translation_md5 = shift;
  throw("Translation_md5 string expected") if (!defined $translation_md5);

  my $translation_stable_id = $self->get_stable_id_for_md5($translation_md5);
  throw("No translation_stable_id for translation_md5 $translation_md5") if (!defined $translation_stable_id);

  my $translation = $self->get_translation($translation_stable_id);
  my $translation_seq = $translation->seq;
  my $transcript = $translation->transcript;
  $self->reverse($transcript->strand < 0);
  my $transcript_stable_id = $transcript->stable_id;

  $self->init_protein_matrix($translation, $translation_md5);

  $self->init_header;

  my $all_triplets = $self->get_triplets($translation_stable_id);

  $self->load_predictions_for_triplets($all_triplets);

  $self->store_protein_matrix($translation_stable_id, $translation_md5) if ($self->{'pipeline_mode'});

  if ($self->{'write_mode'}) {
    if ($translation_seq ne join('', @{$self->amino_acids})) {
      my $fh = FileHandle->new($self->working_dir. "/$translation_stable_id", 'w');
      print $fh "$transcript_stable_id\n$translation_seq\n";
      print $fh join('', @{$self->amino_acids}), "\n";
      $fh->close;
    }
  }
}

=head2 amino_acids
  Arg [1]    : String $aa Amino acid  
  Example    : $runCADDAnnotation->amino_acids('K');
  Description: Holds reference to an array of string. This method is used to collect each
               annotated amino acid. The final array will be compared against the input translation.
  Returntype : Arrayref of string
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub amino_acids {
  my $self = shift;
  my $aa = shift;
  if (defined $aa) {
    push @{$self->{'amino_acids'}}, $aa;
  } else {
    return $self->{'amino_acids'};
  }
}

=head2 annotation_file
  Arg [1]    : String $aa Amino acid  
  Example    : $runCADDAnnotation->amino_acids('K');
  Description: Holds reference to an array of string. This method is used to collect each
               annotated amino acid. The final array will be compared against the input translation.
  Returntype : Arrayref of string
  Exceptions : None
  Caller     : General
  Status     : At Risk
=cut
sub annotation_file {
  my $self = shift;
  return $self->{'annotation_file'};
}

sub annotation_file_version {
  my $self = shift;
  return $self->{'annotation_file_version'};
}

sub assembly {
  my $self = shift;
  return $self->{'assembly'};
}


sub reverse {
  my $self = shift;
  return $self->{'reverse'} = shift if(@_);
  return $self->{'reverse'};
}

sub analysis {
  my $self = shift;
  return $self->{'analysis'} = shift if(@_);
  return $self->{'analysis'};
}

sub working_dir {
  my $self = shift;
  return $self->{'working_dir'};
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
  my $annotation_file = $self->annotation_file;
  open HEAD, "tabix -fh $annotation_file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $header = [split];
  }
  close HEAD;

  $self->header($header);
}

sub parser {
  my $self = shift;
  if (!defined $self->{'parser'}) {
    my $annotation_file = $self->annotation_file;
    $self->{'parser'} = Bio::DB::HTS::Tabix->new(filename => $annotation_file);
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

sub working_dir {
  my $self = shift;
  return $self->{'working_dir'};
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
  my $annotation_file = $self->annotation_file;
  open HEAD, "tabix -fh $annotation_file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $header = [split];
  }
  close HEAD;

  $self->header($header);
}

sub get_stable_id_for_md5 {
  my ($self, $md5) = @_;
  my $var_dba = $self->get_species_adaptor('variation');
  my $get_stable_id_sth = $var_dba->dbc->prepare(qq{
    SELECT  stable_id
    FROM    translation_mapping
    WHERE   md5 = ?
  });

  $get_stable_id_sth->execute($md5);
  my ($stable_id) = $get_stable_id_sth->fetchrow_array;
  return $stable_id if (defined $stable_id);
  return $md5;
}

sub get_translation {
  my $self = shift;
  my $translation_stable_id = shift;
  my $core_type = ($translation_stable_id =~ /^NP|XP/) ? 'otherfeatures' : 'core';
  my $cdba = $self->get_species_adaptor($core_type);
  my $translation_adaptor = $cdba->get_TranslationAdaptor or die "Failed to get translation adaptor";
  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);
  return $translation;
}

sub mutate {
  my $self = shift;
  my $triplet = shift;
  my $reverse = shift;
  my @nucleotides = split('', $triplet);
  my $new_triplets;
  foreach my $i (0 .. $#nucleotides) {
    my $mutations = ['A', 'G', 'C', 'T'];
    $new_triplets = $self->get_mutated_triplets($triplet, $mutations, $i, $new_triplets, $reverse);
  }
  return $new_triplets;
}

sub get_mutated_triplets {
  my $self = shift;
  my $triplet = shift;
  my $mutations = shift;
  my $position = shift;
  my $new_triplets = shift;
  my $reverse = shift;
  foreach my $mutation (@$mutations) {
    my $update_triplet = $triplet;
    if ($reverse) {
      my $reverse_mutation = $mutation;
      reverse_comp(\$reverse_mutation);
      substr($update_triplet, $position, 1, $reverse_mutation);
    } else {
      substr($update_triplet, $position, 1, $mutation);
    }
    $new_triplets->{$triplet}->{$position}->{$mutation} = $update_triplet;
  }
  return $new_triplets;
}

sub get_triplets {
  my $self = shift;
  my $translation_stable_id = shift;
  my $core_type = ($translation_stable_id =~ /^NP|XP/) ? 'otherfeatures' : 'core';
  my $cdba = $self->get_species_adaptor($core_type);
  my $translation_adaptor = $cdba->get_TranslationAdaptor or die "Failed to get translation adaptor";
  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);
  my $slice_adaptor = $cdba->get_SliceAdaptor or die "Failed to get slice adaptor";

  my $transcript = $translation->transcript;
  my $chrom = $transcript->seq_region_name;
  my $start = $transcript->seq_region_start;
  my $end = $transcript->seq_region_end;
  my $strand = $transcript->seq_region_strand;
  my $slice = $slice_adaptor->fetch_by_region('toplevel', $chrom,  $start, $end);
  my $transcript_mapper = $transcript->get_TranscriptMapper();

  my $codonTable = Bio::Tools::CodonTable->new();
  my @all_triplets = ();
  foreach my $i (1 .. $translation->length) {
    my @pep_coordinates = $transcript_mapper->pep2genomic($i, $i);
    my $triplet = '';
    my @coords = ();
    foreach my $coord (@pep_coordinates) {
      my $coord_start = $coord->start;
      my $coord_end = $coord->end;
      next if ($coord_start <= 0);
      my $new_start = $coord_start - $start + 1;
      my $new_end   = $coord_end   - $start + 1;
      my $subseq = $slice->subseq($new_start, $new_end, $strand);
      $triplet .= $subseq;
      push @coords, [$coord_start, $coord_end];
    }
    my $entry = {
      coords => \@coords,
      aa_position => $i,
      chrom => $chrom,
      triplet_seq => $triplet,
    };
    my $aa = $codonTable->translate($triplet);
    if (!$aa) {
      $entry->{aa} = 'X';
    } else {
      $entry->{aa} = $aa;
      my $reverse = ($strand < 0);
      my $new_triplets = $self->mutate($triplet, $reverse);
      $entry->{new_triplets} = $new_triplets;
    }
    push @all_triplets, $entry;
  } 
  return \@all_triplets;
}

sub init_protein_matrix {
  my $self = shift;
  my $translation = shift;
  my $translation_md5 = shift;
  my $pred_matrices = {};
  foreach my $analysis (@{$self->analysis}) {
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
        print $self->working_dir. "/$analysis\_$translation_stable_id", "\n";
        my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, $translation_md5);
        my $debug_data = $self->{debug_data};
        foreach my $i (sort keys %{$debug_data->{$analysis}}) {
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
