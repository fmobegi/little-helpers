#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/..";
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use CoreUtil;
use CoreDB;
use feature qw/switch/;


my ($opt_help, $opt_man, $opt_cmd, $opt_db, $opt_list, $opt_group_id);

GetOptions(
  'help|h'       => \$opt_help,
  'man|m'        => \$opt_man,
  'command|c=s'  => \$opt_cmd,
  'listfile|l=s' => \$opt_list,
  'group_id|g=s' => \$opt_group_id,
  
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit unless ( (scalar @ARGV == 1) and (defined $opt_cmd) );

$opt_db = $ARGV[0];

if (not -e $opt_db) {
    print STDERR "Database does not exist\n";
    return 0;
}
$opt_db = File::Spec->rel2abs($opt_db);

my $util = CoreUtil->new(db => CoreDB->new(db=>$opt_db) );

given ($opt_cmd) {
	when( 'list_groups')           { $util->list_groups }
    when( 'fasta_by_list')         { fasta_by_list($opt_list) }
    when( 'presence_matrix')       { $util->presence_matrix() }
    when( 'fasta_by_group_id')     { fasta_by_group_id($opt_group_id) }
    when( 'protein_length_matrix') { $util->protein_length_matrix }
    default                        { die "unknown command: '$opt_cmd'" }
}

sub fasta_by_list {
	if (not defined $opt_list) {
		pod2usage(-verbose => 1) && exit;
	}
	else {
		$util->fasta_by_list($opt_list);
	}
}

sub fasta_by_group_id {
	if (not defined $opt_group_id) {
		pod2usage(-verbose => 1) && exit;
	}
	else {
		$util->fasta_by_group_id($opt_group_id);
	}
}

__END__

=head1 NAME

=head1 SYNOPSIS

B<pangenomics_util.pl> -command [see command options below] seq.db

 Options:
   -help|h        brief help message
   -listfile|l	  file containing group IDs
   -group_id|g    a group ID
   -command|c     command to execute (see options below)

"B<command>" accepts these options:

=item "list_groups": list groups.

=item "fasta_by_list" (needs -listfile): make fasta files for groups in listfile.

=item "fasta_by_group_id" (needs -group_id): print fasta sequences for a group id.

=item "presence_matrix": print a matrix of presence/absence for each isolate in each group.

	
=head1 DESCRIPTION

=over 8

=item B<-help>

Print a brief help message and exit.

=item B<-man>

Print the manual page and exit.

=back

=head1 VERSION

Version 0.01

=cut
