=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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

=cut


=head1 FinishPhenotypeAnnotation

This module runs at the end of the Phenotype Annotation pipeline and produces 
a summary report to check the results look reasonable.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotationPipeline;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

sub run {
    my $self = shift;
    
    my $hive_dba = $self->dbc;
    
    my $runTime_sth = $hive_dba->prepare(qq[ select timediff(max(when_died), min(when_born)) from worker ]);
    $runTime_sth->execute()||die;
    my $time = $runTime_sth->fetchall_arrayref();

    my $dir =$self->required_param('pipeline_dir');
    open my $report, ">$dir/REPORT_hive_pipe.txt"||die "Failed to open report file for summary info :$!\n";

    print $report "PhenotypeAnnotation pipeline finished! \n";
    print $report "running time: $time->[0]->[0] \n";
    print $report "pipeline_dir: ", $self->required_param('pipeline_dir'), "\n";
    
    close $report;
}

1;

