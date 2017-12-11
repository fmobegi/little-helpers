#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);

use lib "$Bin/../../lib";
use lib "$Bin/../..";

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use CoreStats;
use CoreDB;

my ($opt_help, $opt_type, $opt_db);

GetOptions(
  'help|h'        => \$opt_help,
  'database|d=s'  => \$opt_db,
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 1) && exit unless defined $opt_db;

if (not -e $opt_db) {
    print STDERR "Database does not exist\n";
    return 0;
}
$opt_db = File::Spec->rel2abs($opt_db);

my $coredb = CoreDB->new(db=>$opt_db);


dump_specific_clusters($coredb->get_type_counts_for_clusters);


sub dump_specific_clusters {
    my ($hash_of_cluster) = @_;

         print join "\t", ( 'Analysis',
                            'Peculiarity',
                            'DB_group_id', 
                            'Phenotype1', 
                            'Phenotype1_isolate_count',
                            'Phenotype2', 
                            'Phenotype2_isolate_count',
                            'Gene_product');
         print "\n";

    
    foreach my $cluster_id (keys %{ $hash_of_cluster } ) {
        
        #This array holds one or (usually) two 'Cluster' objects, 
        # one for each phenotype group  ('BAC' or 'MEN') for a given cluster ID.
        #It also holds the counts of members as well as some annotation info.
        my @clusters = @{ $hash_of_cluster->{$cluster_id} };
        my ($cl1, $cl2) = @clusters;
        if ( scalar @clusters == 2 ) { 
            # one the two cluster has 70 members and the total is not 140 (i.e.
            # one is missing some members:
            if ( ($cl1->member_count == 70 or $cl2->member_count == 70) 
                and
                ($cl1->member_count + $cl2->member_count != 140) )
                
            {
                print join "\t", (
                                $cl1->type,
                                'CoreVsNonCore',
                                'GroupID_'.$cluster_id, 
                                $cl1->remark, 
                                $cl1->member_count,
                                $cl2->remark,
                                $cl2->member_count,
                                $cl1->product);
                print "\n";
                
            }
            elsif ( abs($cl1->member_count -  $cl2->member_count) >= 10 ) {
                 print join "\t", (
                                $cl1->type,
                                'CountDiff>10x'
                                ,'GroupID_'.$cluster_id, 
                                $cl1->remark, 
                                $cl1->member_count,
                                $cl2->remark,
                                $cl2->member_count,
                                $cl1->product);
                 print "\n";
            }
        }
        elsif(scalar @clusters == 1) {
            if ($cl1->member_count >=10) { 
                 print join "\t", (
                                $cl1->type,
                                'SingleCluster>10x'
                                ,'GroupID_'.$cluster_id, 
                                $cl1->remark, 
                                $cl1->member_count,
                                0,
                                0,
                                $cl1->product);
                print "\n";
            }
        }
        else {
            die 'Cluster array cannot be empty';
        }
    }
    
}

__END__

=head1 NAME

=head1 SYNOPSIS

B<group_specific.pl> -d seq.db 

 Options:
   -help|h        brief help message
   -database|d    path the to database

=head1 DESCRIPTION

=over 8

=item B<-help>

Print a brief help message and exit.

=back

=head1 VERSION

Version 0.01

=cut
