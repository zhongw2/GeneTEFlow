#!/usr/bin/env perl

# run this script in the directory of the rsem output 
# this script generates mapping stats for rsem output
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


my $sampleinfo = &getsampleinfo( $sampleinfofile );
print join ( "\t", "FastqID", "SampleID", "SampleName", "SampleGroup", "total_read", "mapped_read", "pct_mapped_read" );
print "\n";
foreach my $sample ( sort {$sampleinfo->{$a}->{samplename} cmp $sampleinfo->{$b}->{samplename}} keys %$sampleinfo ) 
{
        my $fastqid = $sample;
        my $sampleid = $sampleinfo->{$sample}->{sampleid};
	my $samplename = $sampleinfo->{$sample}->{samplename};
	my $samplegroup = $sampleinfo->{$sample}->{samplegroup};
        
	my $statfile = "$rsemdir/$sample".".RSEM_Output.cnt";
#	print "$sample\t$statfile\n";
	my $outstat = `head -1 $statfile`;
	my @elms = split ( /\s+/, $outstat );
	my $mapread = $elms[1];
	my $totalread = $elms[3];
#	print "$mapread\t$totalread\n";
	my $pctmap = sprintf( "%.2f", $mapread/$totalread );
	print join ( "\t", $fastqid, $sampleid , $samplename, $samplegroup, $totalread, $mapread, $pctmap );
	print "\n";
}

sub getsampleinfo
{
	my $file = shift;
	my %sampleinfo;
	open ( FILE, $file ) || die " can't open file $file ";
		while ( <FILE> )
		{
			chomp;
			next if $_ =~ /^FastqID/;
			my ( $fastqid, $sampleid, $samplename, $samplegroup, $rest ) = split ( /\t/, $_ );
#                        $fastqid=~ s/$file_extention//g;
			$sampleinfo{$fastqid} = { sampleid=>$sampleid,  samplename=>$samplename, samplegroup=>$samplegroup };
		}
	close FILE;
	return \%sampleinfo;
}
