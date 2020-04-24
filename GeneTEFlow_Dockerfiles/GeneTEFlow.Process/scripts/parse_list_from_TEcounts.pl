#!/usr/bin/perl -w
use strict;


my $usage = "$0 <in_genelist>  <in>\n";

my %hash = ();

my $infile = $ARGV[0];
open(IN, $infile) || die $usage ;
my @array=<IN>;
for(my $j=0; $j<scalar(@array); $j++ )
{
    my $string_array=$array[$j];
    $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
    my $gene_name=$string_array;
    if(!(exists $hash{$gene_name})){$hash{$gene_name}=0;}
    else{print "error\t$gene_name\n";}
 
}

$infile = $ARGV[1];
open(IN, $infile) || die $usage ;

while(<IN>){
    chomp; # avoid \n on last field
    my $line=$_;
    my @value_array=split('\t', $line);
    my $gene_name=$value_array[0];
    my $counts=$value_array[1];
    if(exists $hash{$gene_name}) {$hash{$gene_name}= int($counts);}
}






for(my $j=0; $j<scalar(@array); $j++ )
{
    my $string_array=$array[$j];
    $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
    my $gene_name=$string_array;
    if(exists $hash{$gene_name}){print $gene_name."\t".$hash{$gene_name}."\n";}

}




