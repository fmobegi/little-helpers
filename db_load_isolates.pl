#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use DBI;

my ($opt_help, $opt_man, $opt_remark);
$opt_remark = 'NULL';

GetOptions(
  'help|h'         => \$opt_help,
  'remark|r=s'     => \$opt_remark,
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 1) && exit unless (@ARGV > 1);


my $DB = shift @ARGV;
my $dsn = "DBI:SQLite:$DB";
my $DBH = DBI->connect($dsn, {RaiseError => 1, AutoCommit => 0} );

load_isolates(@ARGV);


sub load_isolates {
    my @isolate_files  = @_;
    
    my $sth = $DBH->prepare("INSERT INTO isolates(sanger_id, remarks) VALUES (?,?)") or die $DBH->errstr;

    my $i = 0;
    my %seen;
    foreach my $filename (@isolate_files) {
        my ($name) = $filename =~ /(.+)\./;
        next if ( $seen{$name} );
        $sth->execute($name, $opt_remark) or die $DBH->errstr;        
        print "Inserting ", ++$i, " rows into database\n";
        $seen{$name} = 1;
    }
    print $i, " rows inserted into database\n";
}



__END__

=head1 NAME

=head1 SYNOPSIS

B<$0> seq.db -r blood_isolates isolate1.fasta isolate2.fasta ...

 Options:
   -help|h        brief help message
   -remark|r      add a 'tag' to the isolates (i.e. blood, CSF, etc...)

=head1 DESCRIPTION

=over 8

=item B<-help>

Print a brief help message and exit.

=back

=head1 VERSION

Version 0.01

=cut
