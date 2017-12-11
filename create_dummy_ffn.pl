#!/usr/bin/env perl
use strict;
use strict;
use Bio::SeqIO;

my $faa = $ARGV[0] or die "No file provided";

my $ffn = $faa;
$ffn =~ s/faa/ffn/;

open my $fh, ">" . $ffn or die $!;

my $seqio = Bio::SeqIO->new(-file => $faa, format=>'fasta');

while (my $seq = $seqio->next_seq) {
    print $fh ">", $seq->id, "\n", "N", "\n";
} 
