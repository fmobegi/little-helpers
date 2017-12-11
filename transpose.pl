#!/usr/bin/perl
##Transpose collumn to rows for huge matrix files, e.g. *.vcf SNP files

use strict ;


sub parseparams {
    my $in  = $ARGV[0] ;
    my $out = $ARGV[1] ;

    if (!-e $in or !defined $in or !defined $out) {
	print "transpose.pl\n" ;
        print "FUNCTION: Transpose a matrix)\nBY: Mobegi_FM (fm3)\n" ;
	print "USAGE: perl transpose.pl inputfile outputfile\n" ;
	exit(2) ;
    }
    return [$in, $out] ;
}

# reads the matrix
sub readtable {
    open(my $handle, $_[0]) or die "could not open file $_[0]" ;
    my @table = () ;
    while (<$handle>) {
	my $line = $_ ;
	chomp($line) ;
	my @spl = split /[\,\t]/, $line ;
	push(@table, \@spl) ;
    }
    close($handle) ;
    return \@table ;
}

# transposes the matrix
sub transpose {
    my $tab  = $_[0] ;
    my @tra  = () ;
    my $mcol = 0 ;	
    my $mcnt = 0 ;
    for(my $cnt=0; $cnt<scalar(@$tab); $cnt++) {
	for(my $col=0; $col<scalar( @{$tab->[$cnt]} ); $col++) {
	    $tra[$col] = [] if (!defined $tra[$col]) ;
	    $tra[$col][$cnt] = $tab->[$cnt][$col] ; 
    	    $mcol = $col if($col > $mcol) ;
	}
	$mcnt = $cnt if($cnt > $mcnt) ;
    }
    return \@tra ;
}

# prints the transposed matrix
sub out {
    my $tra = $_[0] ;
    open(my $handle, ">$_[1]") or die "could not open file $_[1]" ;
    for(my $cnt=0; $cnt<scalar(@$tra); $cnt++) {	
	print $handle $tra->[$cnt][0] ;
	for(my $col=1; $col<scalar( @{$tra->[$cnt]} ); $col++) {
	    if (defined $tra->[$cnt][$col])  {
		print $handle "\t$tra->[$cnt][$col]" ;
	    } else {
		print $handle "\t" ;
	    }
	}
	print $handle "\n" ;
    }
    close($handle) ;
}

sub main {
    my $in  = $_[0] ;
    my $out = $_[1] ;

    my $tab = &readtable($in) ;
    my $tra = &transpose($tab) ;
    &out($tra,$out) ;
}

my $pa = &parseparams() ;
&main($pa->[0], $pa->[1]) ;
