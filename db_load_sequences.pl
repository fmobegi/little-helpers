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
pod2usage(-verbose => 1) && exit unless (@ARGV > 1);


my $DB = shift @ARGV;
my $dsn = "DBI:SQLite:$DB";
my $SCHEMA = Genome::Schema->connect("dbi:SQLite:$DB", { AutoCommit => 0 },);

my %ALLprot;
my %ALLnucl;

load_sequences(@ARGV);

sub load_sequences {

    my @files = @_; #.faa annotation files
    #faa are the protein annotations...
    foreach my $faa ( grep(/faa$/, @files) ) {
        my ($prefix) = $faa =~ /(.+)\./;
        my $ffn = $prefix . '.ffn';
        my $isolate_id = get_isolate_id_by_name($prefix);
        db_insert($faa, $ffn, $isolate_id);
    }
}
print   "Seq:", 'last row'
     , ". Unique_proteins:", scalar keys %ALLprot,
     , ". Unique_dna:", scalar keys %ALLnucl, "\n";

sub db_insert {
    my ($faa, $ffn, $isolate_id) = @_;
    
    my %cur_prot;
    my %cur_prot_annot;
    my %cur_nucl;
   
    $SCHEMA->txn_begin;

    #protein seqs
    my $io = Bio::SeqIO->new(-file => $faa, -format=> 'fasta');
    while(my $seq = $io->next_seq) {
        my $pseq       = $seq->seq;
        my $pid        = $seq->id;
        my $pannot     = $seq->desc;

        $cur_prot{$pid}       = $pseq;
        $cur_prot_annot{$pid} = $pannot;
    }
    
    #nucl seqs
    $io = Bio::SeqIO->new(-file => $ffn, -format=> 'fasta');
    while(my $seq = $io->next_seq) {
        my $nseq       = $seq->seq;
        my $nid        = $seq->id;
        my $nannot     = $seq->desc;
        $cur_nucl{$nid} = $nseq;
    }

    foreach my $id (keys %cur_prot) {
       
        my ($dna_id, $prot_id);
        #have we seen this protein sequence before?
        #if yes, what is the DB id of it?
        if (exists $ALLprot{$cur_prot{$id}}) {
            $prot_id = $ALLprot{$cur_prot{$id}};
        }
        else {
            my $prot_rs = $SCHEMA->resultset('Protein')
                          ->create( { seq => $cur_prot{$id} });
            $prot_id =  $prot_rs->id;
            $ALLprot{$cur_prot{$id}} = $prot_id;
        }
        
        #have we seen this nucleotide sequence before?
        if (exists $ALLnucl{$cur_nucl{$id}} ) {
            $dna_id = $ALLnucl{$cur_nucl{$id}};
        }
        else {
            my $dna_rs = $SCHEMA->resultset('Dna')
                          ->create( { seq => $cur_nucl{$id} });
            $dna_id =  $dna_rs->id;
            $ALLnucl{$cur_nucl{$id}} = $dna_id;
        }

        my $sequence_rs = $SCHEMA->resultset('Sequence')
            ->create( { 
                        dna_id      => $dna_id,
                        protein_id => $prot_id,
                        isolate_id => $isolate_id,
                        product    => $cur_prot_annot{$id},
                      }
                    );
       
       if (not $sequence_rs->id % 1000) {
           print   "Seq:", $sequence_rs->id
                 , ". Unique_proteins:", scalar keys %ALLprot,
                 , ". Unique_dna:", scalar keys %ALLnucl, "\n";
       }
    }

    $SCHEMA->txn_commit;
}

sub get_isolate_id_by_name {

    my $name = shift;
    my $isolate_id;    
    eval {
        $isolate_id = $SCHEMA->resultset('Isolate')->search({ sanger_id => $name })->first->id;
    };
    if ($@) {
        die "Looks like isolate with name $name is not in our DB, dying!". $@;
    }

    return $isolate_id;
}

__END__

=head1 NAME

=head1 SYNOPSIS

B<$0> [SQLite DB] *.faa ...

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
