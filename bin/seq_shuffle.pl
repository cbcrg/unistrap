#!/usr/bin/perl
use List::Util 'shuffle';

($aln_file, $out_file, $seed, $replicate_num)=@ARGV;

#_________ Read MSA file and create column hashes _________#


%ALNhash=(); 
%NAMhash=();
$refNum=0;
$count=0;
$/=">";
open ALNseq, $aln_file;

$count=0;
while (<ALNseq>)
{   
    $sequence=""; $title="";
    $entry=$_;
    chop $entry;
    $entry= ">"."$entry";
    $entry=~/>(.+?)\n(\C*)/g;
    $title=$1;$sequence=$2;
    $titleNoU=$title; #$titleNoU=~s/U//g;
    $sequence=~s/\n//g; #$sequence=~s/-//g;
   if ($titleNoU ne "" && $sequence ne "")
   {     $count++;
	 $ALNhash{$title}=$sequence;
 	 $NAMhash{$count}=$title;
	 #print  "$count : ".$title."  $sequence\n";
   }
} 

$/="\n"; #print "Reading of MSA file DONE!!\n"; 

$seed_to_use=$replicate_num+$seed;

for($i=0; $i<$replicate_num; $i++)
{
    srand($seed_to_use+$i); 
 
    $num=$i+1;
    open OUT, ">", "${out_file}_$num.fa" or die $!; 

    @values = shuffle(values %NAMhash);
    $NAMhash{$_} = shift @values for keys %NAMhash;

    foreach ( keys %NAMhash) {
      print OUT ">$NAMhash{$_}\n$ALNhash{$NAMhash{$_}}\n";
    }
    close(OUT);
}


