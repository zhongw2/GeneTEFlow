#!/usr/bin/env perl

# this script generates mapping stats for STAR mapping

use strict;
use warnings;

if ( @ARGV < 1 )
{
	print "Usage: <rsemdir> <sampleinfo file> > <all.sample.map.stat.out> ";
	exit(1);
}

my $rsemdir = $ARGV[0];
my $sampleinfofile = $ARGV[1];
my $file_extention = $ARGV[2];


#my $sampleinfo = &getsampleinfo( $sampleinfofile );
my ($sampleinfo,@sampleinfo) = &getsampleinfo( $sampleinfofile );
print join ( "\t", "FastqID", "SampleID", "SampleName", "SampleGroup", "Total_read","Unique_mapped_read","Multiple_mapped_read", "Total_mapped_read", "Pct_mapped_read" );
print "\n";

#foreach my $sample ( sort {$sampleinfo->{$a}->{samplename} cmp $sampleinfo->{$b}->{samplename}} keys %$sampleinfo ) 
foreach my $sample (@sampleinfo)
{
        my $fastqid = $sample;
        my $sampleid = $sampleinfo->{$sample}->{sampleid};
	my $samplename = $sampleinfo->{$sample}->{samplename};
	my $samplegroup = $sampleinfo->{$sample}->{samplegroup};
        
	my $statfile = "$rsemdir/$sample".".RSEM_Output.STAR.log";
#	print "$sample\t$statfile\n";
#	my $outstat = `head -1 $statfile`;
#	my @elms = split ( /\s+/, $outstat );
#	my $mapread = $elms[1];
#	my $totalread = $elms[3];
#	print "$mapread\t$totalread\n";
        my $totalread = `cat $statfile | grep "Number of input reads"|cut -f 2`;
        my $mapread_uniq = `cat $statfile|grep "Uniquely mapped reads number"|cut -f 2`;
        my $mapread_multi = `cat $statfile|grep "Number of reads mapped to multiple loci"|cut -f 2`;
        chomp $totalread;
        chomp $mapread_uniq;
        chomp $mapread_multi;
        my $mapread = $mapread_uniq + $mapread_multi;
  	my $pctmap = sprintf( "%.2f", $mapread/$totalread );
	print join ( "\t", $fastqid, $sampleid, $samplename, $samplegroup, $totalread, $mapread_uniq, $mapread_multi, $mapread, $pctmap );
	print "\n";
}

sub getsampleinfo
{
	my $file = shift;
	my %sampleinfo;
        my @sampleinfo;
	open ( FILE, $file ) || die " can't open file $file ";
		while ( <FILE> )
		{
			chomp;
			next if $_ =~ /^FastqID/;
			my ( $fastqid,$sampleid, $samplename, $samplegroup, $rest ) = split ( /\t/, $_ );
#                        $fastqid=~ s/$file_extention//g;
			$sampleinfo{$fastqid} = { sampleid=>$sampleid, samplename=>$samplename, samplegroup=>$samplegroup };
                        push @sampleinfo, $fastqid;
		}
	close FILE;
#	return \%sampleinfo;
        return (\%sampleinfo , @sampleinfo)  ;
}
