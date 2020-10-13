#!/usr/bin/env perl
# this scripts generates count matrix input for deseq analysis or TPM matrix
use strict;
use warnings;

my $rsemdatadir=$ARGV[0];  # top dir directory
my $sampleinfofile=$ARGV[1];
my $datatype=$ARGV[2];
my $offsite = -1; # for new file formats; expected_count or TPM column
 

if ( @ARGV < 3 )
{
	print "Usage $0 <rsem results dir> <sampleinfofile> <type of data to retrive (counts, TPM) \n";
	exit(1);
}

my ($sampleinfo,@sampleinfo) = &getsampleinfo( $sampleinfofile );
#my @sampleinfo = &getsampleinfo( $sampleinfofile );
my @generesultfile;
#foreach my $sample ( sort {$sampleinfo->{$a}->{samplename} <=> $sampleinfo->{$b}->{samplename}} keys %$sampleinfo ) 
#foreach my $sample ( keys %$sampleinfo ) 
#foreach my $sample (@sampleinfo)
#foreach my $sample ( sort {$sampleinfo->{$a}->{samplename} cmp $sampleinfo->{$b}->{samplename}} keys %$sampleinfo ) 
foreach my $sample (@sampleinfo)
{
        my $fastqid = $sample;
        my $sampleid = $sampleinfo->{$sample}->{sampleid};
	my $samplename = $sampleinfo->{$sample}->{samplename};
	my $samplegroup = $sampleinfo->{$sample}->{samplegroup};
        

	my $generesult = "$rsemdatadir/$sample".".TEcount.txt";
	
	my $degeneresultln = "$samplename.TE.results";
#	print "$generesult\t$degeneresultln\n";
	unless ( -e  $degeneresultln )
	{
		`ln -s $generesult $degeneresultln`;
	}
#        my $degeneresultln = $generesult;
	push @generesultfile,  $degeneresultln;
}
my $allresultfile = join ( " " , @generesultfile );
if ( $datatype eq "counts" )
{
	my $outputfile = "all.sample.Counts.TE.results.txt";
#	`~/rsem-generate-data-matrix $allresultfile  > $outputfile  ` ; 

        $offsite = 1;       # for new file formats, expected_count      
  
        rsemGenerateDataMatrix($allresultfile,$outputfile)
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
                        my ( $fastqid, $sampleid, $samplename, $samplegroup, $rest ) = split ( /\t/, $_ );
                        $sampleinfo{$fastqid} = { sampleid=>$sampleid, samplename=>$samplename, samplegroup=>$samplegroup };
						push @sampleinfo, $fastqid;
                }
        close FILE;
#        return @sampleinfo;
        return (\%sampleinfo , @sampleinfo)  ;
}


sub rsemGenerateDataMatrix
{
    my ($allresultfile_string,$outputfile) = @_; 
    my @allresultfileArray=split(' ',$allresultfile_string);

#    my $offsite = 5; # for new file formats; TPM column

    my $line;
    my $n = scalar(@allresultfileArray);
    my $M = -1;
    my @matrix = ();

    my @ids_arr = ();

    for (my $i = 0; $i < $n; $i++) {
         my (@ids, @ecs) = ();
         &loadData($allresultfileArray[$i], \@ecs, \@ids);

        if ($M < 0) { 
            $M = scalar(@ids); 
            @ids_arr = @ids;
        }
        elsif (!&check($M, \@ids_arr, \@ids)) { 
            print STDERR "Number of lines among samples are not equal!\n"; 
            exit(-1); 
        }

        
        my @filename_array=split(".genes.results",$allresultfileArray[$i]);
        my $real_filename = $filename_array[0];
         
        @ecs = ($real_filename, @ecs);
#        @ecs = ($allresultfileArray[$i], @ecs);
        push(@matrix, \@ecs);
    }


    @ids_arr = ("", @ids_arr);
    @matrix = (\@ids_arr, @matrix);

    open(my $fh, '>', $outputfile);
    for (my $i = 0; $i <= $M; $i++) {
#         print $fh "$allresultfileArray[$i]\n";
         for (my $j = 0; $j < $n; $j++) { 
              print $fh "$matrix[$j][$i]\t"; 
         }
         print $fh "$matrix[$n][$i]\n";
    }

    close $fh;
}


# 0, file_name; 1, reference of expected count array; 2, reference of transcript_id/gene_id array
sub loadData {

    #my $offsite = 5; # for new file formats; TPM column

    if($offsite < 0) {
            print STDERR "Double check the value of offsite!\n";
            exit(-1);
    }


    open(INPUT, $_[0]);
    my $line = <INPUT>; # The first line contains only column names
    while ($line = <INPUT>) {
        chomp($line);
        my @fields = split(/\t/, $line);
        push(@{$_[2]}, $fields[0]);
        push(@{$_[1]}, $fields[$offsite]);
    }
    close(INPUT);

    if (scalar(@{$_[1]}) == 0) {
        print STDERR "Nothing is detected! $_[0] may not exist or is empty.\n";
        exit(-1);
    }
}


#0, M; 1, reference of @ids_arr; 2, reference of @ids
sub check {
    my $size = $_[0];
    for (my $i = 0; $i < $size; $i++) { 
        if ($_[1]->[$i] ne $_[2]->[$i]) {
            return 0;
        }
    }
    return 1;
}


