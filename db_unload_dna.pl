#!/usr/bin/env perl
use strict;
use warnings;
use lib '..';
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Genome::Schema;

my ($opt_help, $opt_man);

GetOptions(
  'help|h'         => \$opt_help,
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 1) && exit unless (@ARGV == 1);


my $DB = shift @ARGV;
my $dsn = "DBI:SQLite:$DB";
my $SCHEMA = Genome::Schema->connect("dbi:SQLite:$DB", { AutoCommit => 0 },);

my $rs1 = $SCHEMA->resultset('Dna');

open my $fh1, ">dna.fasta" or die $!;
while(my $seq = $rs1->next) {
    my $st = $seq->seq;
    $st =~ s/\*$//;
    print $fh1 ">d", $seq->id, "\n", $st, "\n";
}
        
__END__

=head1 NAME

=head1 SYNOPSIS

B<$0> [SQLite DB] 

 Options:
   -help|h        brief help message

=head1 DESCRIPTION

=over 8

=item B<-help>

Print a brief help message and exit.

=back

=head1 VERSION

Version 0.01

=cut
