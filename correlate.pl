#!/usr/bin/perl -w
use strict;
use Math::Libm ':all';

## F.M. Mobegi, PhD
## 03/05/2017
## Reference is the first line in the file. 
my $usage = "Usage: $0 pearson/spearman rows/cols infile\n";

if (scalar (@ARGV) != 3) {
	die $usage; }
if (($ARGV[0] ne "pearson") && ($ARGV[0] ne "spearman")) {
	die $usage; }
if (($ARGV[1] ne "rows") && ($ARGV[1] ne "cols")) {
	die $usage; }
if (! (open (I, "<$ARGV[2]"))) {
	die "Can't open $ARGV[2]: $!\n"; }

my %data = ();
my $line = <I>;
chomp ($line);
$line =~ s/^.*?\t/\t/;
if ($line =~ /^[\t\-\d\.eE\+]+\Z/) {
	print STDERR "Warning: assuming these are column heads:\n$line\n"; }
my @colheads = split /\t/o, $line;
if (scalar @colheads == 1) {
	die "Columns should be separated by tabs\n$usage"; }
while ($line = <I>) {
	chomp ($line);
	my @a = split /\t/o, $line;
	if ($a[0] =~ /^[\-\d\.eE\+]+\Z/o) {
		print STDERR "Warning: assuming this is a row head: $a[0]\n"; }
	for (my $col = 1; $col < scalar @a; ++$col) {
		if ($ARGV[1] eq "cols") {
			push (@{$data{$colheads[$col]}}, $a[$col]); }
		else {
			push (@{$data{$a[0]}}, $a[$col]); } } }

foreach my $id1 (sort keys %data) {
	ID2: foreach my $id2 (sort keys %data) {
		if ($id2 le $id1) {
			next ID2; }
		print STDOUT "$id1\t$id2\t";
		my $correlation = 0;
		if ($ARGV[0] eq "pearson") {
			$correlation = &pearson(\@{$data{$id1}}, \@{$data{$id2}}); }
		else {
			$correlation = &spearman(\@{$data{$id1}}, \@{$data{$id2}}); }
		print STDOUT "@{$correlation}[0]\t@{$correlation}[1]\n"; } }

sub pearson () {
	my @array1 = @{$_[0]};
	my @array2 = @{$_[1]};
	
	my $n = scalar @array1;
	if ($n != scalar @array2) {
		die "ERROR: arrays passed to subroutine 'pearson' have different lengths\n"; }
	if ($n == 0) {
		die "ERROR: arrays passed to subroutine 'pearson' have zero length\n"; }

	my $Sxy = 0;
	my $Sx = 0;
	my $Sy = 0;
	my $Sxx = 0;
	my $Syy = 0;
	my %values = ();
	for (my $i = 0; $i < $n; ++$i) {
		++$values{"1"}{$array1[$i]};
		++$values{"2"}{$array2[$i]};
		$Sxy += $array1[$i] * $array2[$i];
		$Sx += $array1[$i];
		$Sy += $array2[$i];
		$Sxx += $array1[$i] ** 2;
		$Syy += $array2[$i] ** 2; }
	
	if ((scalar (keys %{$values{"1"}}) == 1) || (scalar (keys %{$values{"2"}}) == 1)) {
		print STDOUT "ERROR: one of the arrays passed to subroutine 'pearson' contains completely monotonous values\n";
		die; }
	
	my @output = ();
	if (sqrt (($Sxx - (($Sx * $Sx) / $n)) * ($Syy - (($Sy * $Sy) / $n)))) {
		$output[0] = ($Sxy - (($Sx * $Sy) / $n)) / (sqrt (($Sxx - (($Sx * $Sx) / $n)) * ($Syy - (($Sy * $Sy) / $n))));
		$output[1] = 1 - erfc (abs ($output[0]) * sqrt ($n) / sqrt (2)); }
	return (\@output); }

sub spearman () {
	my @array1 = @{$_[0]};
	my @array2 = @{$_[1]};

	my $n = scalar @array1;
	if ($n != scalar @array2) {
		die "ERROR: arrays passed to subroutine 'spearman' have different lengths\n"; }
	if ($n == 0) {
		die "ERROR: arrays passed to subroutine 'spearman' have zero length\n"; }

	my @sort1 = sort { $a <=> $b; } @array1;
	my %ranks1 = ();
	for (my $i = 0; $i < $n; ++$i) {
		my $rank = $i;
		while (($i < $n - 1) && ($sort1[$i + 1] == $sort1[$i])) {
			++$i; }
		for (my $j = $rank; $j <= $i; ++$j) {
			$ranks1{$sort1[$j]} = (($rank + $i) / 2) + 1; } }
	my @sort2 = sort { $a <=> $b; } @array2;
	my %ranks2 = ();
	for (my $i = 0; $i < $n; ++$i) {
		my $rank = $i;
		while (($i < $n - 1) && ($sort2[$i + 1] == $sort2[$i])) {
			++$i; }
		for (my $j = $rank; $j <= $i; ++$j) {
			$ranks2{$sort2[$j]} = (($rank + $i) / 2) + 1; } }
	my $Sx_y2 = 0;
	for (my $i = 0; $i < $n; ++$i) {
		$Sx_y2 += ( $ranks1{$array1[$i]} - $ranks2{$array2[$i]} ) ** 2; }
	my @output = ();
	$output[0] = 1 - ((6 * $Sx_y2) / (($n ** 3) - $n));
	$output[1] = "P-value not implemented for Spearman correlation";
	return (\@output); }

sub mean () {
	my @array = @{$_[0]};
	my $sum = 0;
	foreach my $i (@array) {
		$sum += $i; }
	return ($sum / (scalar @array)); }

