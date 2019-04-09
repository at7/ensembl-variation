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

package Bio::EnsEMBL::Variation::Utils::RunCADDAnnotationUtils;

use Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotationUtils;
our @ISA = ('Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotationUtils');

my $CADD_CUTOFF = 0.5;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  if (! grep {$_ eq $self->annotation_file_version} ('v1.3', 'v1.4', 'v1.5')) {
    die "CADD version " . $self->annotation_file_version . " is not supported.";
  }

  $self->analysis([qw/cadd/]);

  return $self;
}


sub load_predictions_for_triplets {
  my $self = shift;
  my $all_triplets = shift; 
  foreach my $entry (@$all_triplets) {
    my $aa = $entry->{aa};
    $self->amino_acids($aa);
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
        my $cadd_phred = $data->{'PHRED'};
        next if ($alt eq $ref);
        my $nucleotide_position = ($self->reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $mutated_triplet =  $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $self->codon_table->translate($mutated_triplet);
        $self->add_predictions($data, $i, $mutated_aa);
      }
    }
  }
} 

sub add_predictions {
  my ($self, $data, $i, $mutated_aa) = @_;
  if ($data->{'PHRED'} ne '.') {
    my $prediction = ($data->{'PHRED'} >= $CADD_CUTOFF) ? 'likely deleterious' : 'likely benign';
    $self->add_prediction($i, $mutated_aa, 'cadd', $data->{'PHRED'}, $prediction);
  }
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

1;
