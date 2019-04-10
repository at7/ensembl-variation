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

This module contains helper methods for managing database connections
and API adaptors. 

=cut


package Bio::EnsEMBL::Variation::Utils::BaseDatabaseUtils;

use strict;
use warnings;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Registry;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($species, $registry_file) = rearrange([qw(SPECIES REGISTRY_FILE)], @_);
  my $self = bless {
      'species' => $species,
      'registry_file' => $registry_file,
  }, $class;

  return $self;
}

sub adaptor {
  my $self = shift;
  return $self->{'adaptor'};
}

sub species {
  my $self = shift;
  return $self->{'species'};
}

sub registry_file {
  my $self = shift;
  return $self->{'registry_file'};
}

sub get_species_adaptor {
  my ($self, $group) = @_;
  my $species = $self->species;
  return $self->get_adaptor($species, $group);
}

sub get_adaptor {
  my ($self, $species, $group) = @_;
  my $dba;
  eval {
    $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
  };

  unless (defined $dba) {
    $self->_load_registry();
    $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
  }

  unless (defined $dba) {
    die "Failed to a get DBA for $species and group $group";
  }
  return $dba;
}

sub _load_registry {
  my ($self) = @_;
  my $reg_file = $self->registry_file;
  die("ERROR: Registry file $reg_file not found\n") unless -e $reg_file;
  Bio::EnsEMBL::Registry->load_all($reg_file, 0, 1);
  return;
}


1;