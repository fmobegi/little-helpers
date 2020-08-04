#!/usr/bin/perl -w

my %hash=();
my %id=();
my @id_array=();
my $option=0;
my %Tag=();
my %Tagc=();
my %g_id=();
#gene55970: 880/1071	Percentage:82.1662	Depth:4
#gene1537538: 0/1011	Percentage:0	Depth:nan
#gene55971: 0/708	Percentage:0	Depth:nan
my $List = $ARGV[0];
open ($fh,'<', $List) or die "Could not open file '$List' $!";
open(my $out, '>', 'Merged_Data.txt') or die "Could not open file 'Merged_Data.txt' $!";
print $out "gene\tlength\t";
while(my $line = <$fh>){
   chomp ($line);
      my @DATA=split(/\s+/,$line);
print "Check the File format: Two columns separated by Tab\n\n" if scalar(@DATA)<2;
exit if scalar(@DATA)<2;
	if ($DATA[1] =~ /^ID$/){
		unshift @id_array,$DATA[1];
		}
	else{
   push @id_array,$DATA[1];
		}
#   print $sample_id,"\t";
   open (FILE,$DATA[0]) or die $!;
   while(<FILE>){
     chomp;
     my @data = split(/\t/,$_);
	if ($data[0] =~ /\.hmm$|\_/){
		$data[2] =~ s/\_\d+//;
		$Tagc{$data[2]} = $data[0];# rwmove "$data[0]" and replace with "$_" if wanted all description columns
		$option =1;
		next;
	}
my ($id,$length) = split(/\:\s+/,$data[0]);
$length =~ s/\d+\///;
$data[2] =~ s/Depth://;
$data[2] =~ s/nan/0/;
	$Tag{$id}=1;
     $g_id{$id}=$length;
     $hash{$DATA[1]}{$id}=$data[2];

     }
   }

%Tag = %Tagc if $option==1;
print $out "$_\t" foreach @id_array;
print $out "\n";
shift @id_array if $option==1;
foreach $gene_id (sort keys%Tag){
if ($Tag{$gene_id} =~ /^1$/){
print $out $gene_id,"\t",$g_id{$gene_id},"\t";
}
else{
  print $out $gene_id,"\t",$g_id{$gene_id},"\t",$Tag{$gene_id},"\t";
}
  foreach my $line (@id_array){
    if (exists $hash{$line}{$gene_id}){
      print $out $hash{$line}{$gene_id},"\t";
      }
    else{
      print $out "0\t";
      }
    }
    print $out "\n";
    }
    