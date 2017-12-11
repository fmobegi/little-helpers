#!/usr/bin/env perl

use strict;
use warnings;
use IPC::Cmd qw(run);
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

my ($opt_help, $opt_db, $opt_stage, $opt_blast);
my ($opt_mcl_i, $opt_meta) = (4, 'N.A.');

GetOptions(
  'help|h'        => \$opt_help,
  'stage|s=s'     => \$opt_stage,
  'database|d=s'  => \$opt_db,
  'metainfo|m=s'  => \$opt_meta,
  'blastfile|b=s' => \$opt_blast,
  'inflation|i=i' => \$opt_mcl_i,
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 1) && exit unless args_are_okay();


$opt_stage =~ /^db/i and run_database();
$opt_stage =~ /^load/i and run_load();
$opt_stage =~ /^pre/i  and run_pre_blast();
$opt_stage =~ /^post/i and run_post_blast();

sub run_pre_blast {
    run_this("mkdir set_prot");
    chdir 'set_prot';
    run_this("db_unload_protein.pl $opt_db");
    run_this('makeblastdb -in protein.fasta -input_type fasta -parse_seqids -dbtype prot');
    chdir "..";
    mkdir 'set_dna';
    chdir 'set_dna';
    run_this("db_unload_dna.pl $opt_db");
    run_this('makeblastdb -in dna.fasta -input_type fasta -parse_seqids -dbtype nucl');
}

sub run_post_blast {

    run_this("cut -f 1,2,11 $opt_blast > seq.abc");
    run_this("mcxload -abc seq.abc -write-tab seq.dict -o seq.mci --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)'");
    run_this("mcl seq.mci -I $opt_mcl_i");
    run_this('mcxdump -icl out.seq.mci.I40 -o dump.seq.mci.I140 -tabr seq.dict');
    run_this("db_load_groups.pl $opt_db dump.seq.mci.I140");
}

sub run_load {
    run_this("db_load_isolates.pl -remark $opt_meta $opt_db @ARGV");
    run_this("db_load_sequences.pl $opt_db @ARGV");
}

sub run_database {
    run_this("sqlite3 $opt_db < /nfs/users/nfs_f/fy2/software/CoreGenomeUnofficial/config/core.sql");
}

sub args_are_okay {

    return 0 unless (defined $opt_stage and defined $opt_db);
    return 0 unless $opt_stage =~ /^(db|load|pre|post)/i;
    
    if ($opt_stage !~ /^db/) {
        if (not -e $opt_db) {
            print STDERR "Database does not exist\n";
            return 0;
        }
    } 
    else {
        if (-e $opt_db) {
            print STDERR "Database already exist\n";
            return 0;
        }
    }

    if ($opt_stage =~ /^post/i) {
        if (not defined $opt_blast) {
            print STDERR "You need to provide a BLAST tabular file\n";
            return 0;
        }
    }

    if ($opt_stage =~ /^load/i) {
        if (@ARGV < 2) {
            print STDERR 'You need to provide .ffn and .faa files.';
            return 0;
        }
    }

    $opt_db = File::Spec->rel2abs($opt_db);
    $opt_blast = File::Spec->rel2abs($opt_blast);
    return 1;
}

sub run_this {
    my $c = shift;
    my $t  = scalar localtime;
    #can_run($c) or print STDERR "ERROR: cannot run $c\n $!" and exit 1;
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
        run( command => $c, verbose => 1 );
    if ( $success ) {
        print "Ran at $t\n";
    }
    else {
        print  "ERROR: $c at $t\n"
              . Dumper $error_message
              . Dumper $stderr_buf;
        exit 1;
    }
}

__END__

=head1 NAME

=head1 SYNOPSIS

B<run_core.pl> -d seq.db -t [db|dna|protein] -s [db|load|preblast|postblast] 

 Options:
   -help|h        brief help message
   -stage|s       REQUIRED [db|load|preblast|postblast]
   -database|d    REQUIRED database (A new DB name if stage = 'db')
   -inflation|i   Default = 4, mcl inflation value
   -metainfo|m    OPTIONAL, use for 'load' stage
   -blastfile|b   REQUIRED for 'postblast'stage

=head1 DESCRIPTION

=over 8

=item B<-help>

Print a brief help message and exit.

=back

=head1 VERSION

Version 0.01

=cut
