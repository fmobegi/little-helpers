#!/usr/bin/perl -w

### This script creates multiple submissin scripts for tango ERSA servers and submits them to sbatch queue.

###INPUT ARGUMENTS ######################################################
my $inputdir = "input_directory"; 
my $outdir = "output_directory";

### CHANGE VALUES AS PER NEED ##########################################

my $jobname = "name";
my $ppn = 0-64; 
my $mem = "xxxgb"; 
my $walltime = "100:00:00"; 
my $modules = "
module load moduleName
module load moduleName
"; 

my $cmdline = "
kraken --preload --threads 16 --fasta-input --output outdir\/tag --db /data/genomicsdb/minikraken/ input1
";

$inputdir =~ s/\/$//; 
$outdir =~ s/\/$//; 
my $subdir = "$inputdir\/subscripts"; 
my $cmd = $cmdline;

if (-e $subdir and -d $subdir) {
    print "$subdir exists :)\n";
    }
else {
    mkdir $subdir;
    }

if (-e $outdir and -d $outdir) {
    print "$outdir exists :)\n";
    } 
else {
    mkdir $outdir;
    }

opendir(DIR,$inputdir) or die 
print "Provide Folder Name conatning pair-end barcode and sequencing files\n$!"; 

my %files=(); 
my %file_tag=(); 

while(my $file = readdir(DIR)){
    #	print $file,"\n";
    if ($file =~ m/(\w+)\.(\w+)$/){
    #	$tag{$1}++; $file_tag{$1}{name}{$file}=1;
        my $tag = $1; 
        $tag =~ s/_[L|R].*//; 
        $file_tag{$tag}++; 
        $files{$tag}{name}{$file}=1;
        }
    }

my $c=0; 
my @files_sub=(); 

foreach my $t1(keys%file_tag){ 
    my @file_name=(); 
    $c++; 
    push 
    @file_name,"$inputdir\/$_" foreach (keys%{$files{$t1}{name}});
    #print $t1,"\t",$file_name[0],"\t",$file_name[1],"\n";
    $cmdline =~ s/outdir/$outdir/; 
    $cmdline =~ s/tag/$t1/; 
    $cmdline =~ s/input1/$file_name[0]/; 
    $cmdline =~ s/input2/$file_name[1]/; 

    my $line = "\#!/bin/csh

    ### Job Name
    #SBATCH --job-name=$jobname-$c

    ### Set email type for job
    ### Accepted options: NONE, BEGIN, END, FAIL, ALL
    #SBATCH --mail-type=ALL

    ### email address for user
    #SBATCH --mail-user=ll\@sahmri.com

    ### Queue name that job is submitted to
    #SBATCH --partition=tango

    ### Request nodes
    #SBATCH --ntasks=$ppn
    #SBATCH --mem=$mem
    #SBATCH --time=$walltime

    echo Running on host `hostname`
    echo Time is `date`

    #module(s) if required module load application_module

    $modules\n
    $cmdline \n"; 

    open (SUBFILE,">$subdir\/$jobname\_$c\.sub"); 
    print SUBFILE $line,"\n"; 
    push @files_sub,"$subdir\/$jobname\_$c\.sub";
    $cmdline= $cmd;
}
#exit;

foreach $n (0 ..scalar(@files_sub)-1){
    my $sub_script = $files_sub[$n];
    my $jobname = $sub_script;
    $jobname =~ s/$subdir\///;
    $jobname =~ s/\.sub//;
    $jobname =~ s/\_/-/;
    `sbatch -e $jobname\.o%J $sub_script`;
    ##`sbatch -C "dc:ep" -e $jobname\.o%J $sub_script`; ## force jobs to current queue (sbatch partition) only
}

