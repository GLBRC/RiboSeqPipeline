#!/usr/bin/perl

use strict;

# This function loads gene vectors from Mike Place analysis of Ribo-seq counts
# His file is space delim, rows we care about are Y\d\d\n\n\n\d.+ followed by vector of
# read counts from -72 bp of the start codon to +60 beyond the end of the gene.  

# This script will normalize each vector according to the method of Stein et al. (Frydman lab):
# discount first and last 20 codons of each gene and take the mean from 20*3 = 60 nt into the gene and
# -60 nt from the end of the gene.
#
# Mike's script gives count data from -72 to +60 of the gene ... thus, to calculate
# mean across gene body for normalization:  sum from position 132 of the vector to (end - 120) of the vector, 
# then take average across the measured codons.
# Use this script to write out:  the sum, length of the remaining vector (so can calculate avg
# reads per bp to filter as Stein did), and avg reads per gene-body ... then write out data as % reads PER CODON normalized to
# center-gene body mean reads

my ($vectorfile, $savefile) = @ARGV;
if (scalar @ARGV < 2)  {
    print "usage:  Gene-vector file (-72 to +60 each ORF), savefile \n";
    die;
}

#vectorfile format:  [0] = UID, all other positions are counts from -72 to +60 of each gene
my @Vectors = load_SpaceDelimFile("$vectorfile");

#savefile format:  write out two:  SF is just bp-normalized counts from -72 to + 60
# SF2 is normlized CODON counts from 0 (ATG, = positions -13, -12, -11) to end of the gene

my $savefile1 = $savefile . "_Norm_bpCounts.txt";
open (SF, ">", $savefile1 || die "cannot open $savefile1\n");
print SF "UID\tSumGene-bodyCounts\tPositionsPerGeneBody\tAvgCounts\t";
my $v = -72;  # write out vector header
for (my $r=0; $r<4000; $r++)  {  # hard code header length for simplicity
    print SF "$v\t";
    $v++; 
}
print SF "\n";
my $savefile2 = $savefile . "_Norm_CodonCounts.txt";
open (SF2, ">", $savefile2 || die "cannot open $savefile2\n");
print SF2 "UID\tSumGene-bodyCounts\tCountsPerGeneBody\tAvgCounts\t";
my $v = 0;  # write out vector header for codon position
for (my $r=0; $r<2000; $r++)  {  # hard code header length for simplicity
    print SF2 "$v\t";
    $v++; 
}
print SF2 "\n";

#  Add SF 3:  add a pseudocount to each codon position, still normalize to mean of gene body
my $savefile3 = $savefile . "_Pseudcount-Norm_CodonCounts.txt";
open (SF3, ">", $savefile3 || die "cannot open $savefile3\n");
print SF3 "UID\tSumGene-bodyCounts\tCountsPerGeneBody\tAvgCounts\t";
my $v = 0;  # write out vector header for codon position
for (my $r=0; $r<2000; $r++)  {   # hard code header length for simplicity
    print SF3 "$v\t";
    $v++; 
}
print SF3 "\n";

# now for every gene in the file, calculate mean reads in gene body as described above;
# will write that out and also use it to calcullate normalized reads

for (my $r=0; $r<scalar(@Vectors); $r++)   {  # forevery row on the gene file

    my $gene = $Vectors[$r][0];
    my $genebodystart = (72+60);  # vector starts at -72 and then want to start 60 nt into gene
    print SF "$gene\t";
    print SF2 "$gene\t";
    print SF3 "$gene\t";
    
    my ($sum, $pseudosum, $counts);  # for calculating average, also for pseudocount output

        for (my $c = $genebodystart; $c< ( scalar(@{ $Vectors[$r]} ) - 120);  $c++)   {
            $sum += $Vectors[$r][$c];
            $counts++;  # number basepairs per gene body
        }
        my $avg;  #average read counts per gene body
        if ($counts>0)  {
                $avg = $sum / $counts;
        }

        print SF "$sum\t$counts\t$avg\t";
        print SF2 "$sum\t$counts\t$avg\t";
        print SF3 "$sum\t$counts\t$avg\t";

        # adjust sum for pseudo counts:  below add 1 pseudo count per codon, thus here
        # pseudosum = sum + (gene body length /3) = the # of codons in the gene

        my $genebodlength = ( scalar(@{ $Vectors[$r]} ) - 120) - $genebodystart;  # length of central region
        $pseudosum = $sum + $genebodlength;
        #print "pseudosum $pseudosum = $sum + $genebodlength\n";

        # write out norm counts per position first in original vector format
        # easier to loop in two separate stages even if inelegant
        if ($avg > 0)  {
            for (my $c = 1; $c<scalar(@{ $Vectors[$r] }); $c++)  { #for the length of the vector
                my $val = $Vectors[$r][$c];
                $val = $val / $avg;

                print SF "$val\t";
            }

            print SF "\n";

            # now by codon:  start at vector position -12 (position 60) and sum counts over every 3 positions, this is the ATG
            for (my $cc = 61;  $cc< scalar(@ { $Vectors[$r] }); $cc++)   { # ** but gene name is cell [0] so +1
                my $codonsum;
                      
                #print "starting with $cc position\n";
             
                for (my $codon = $cc; $codon<($cc+ 3); $codon++)   {
                    $codonsum += $Vectors[$r][$codon];
                }
                $cc += 2;  #  bump up $c by two hopefully, then another bump up will start the next codon
                my $numcodons = $counts / 3;
                my $codonavg = $sum / $numcodons;
                my $codonval = $codonsum / $codonavg;  # now norm to avg # reads per codon in the gene body
                my $pseudoavg = $pseudosum / $numcodons;  
                my $pseudoval = ($codonsum + 1) / $pseudoavg;  # now norm to pseudo-count average not sum

                print SF2 "$codonval\t";
                print SF3 "$pseudoval\t";
            }
            print SF2 "\n";
            print SF3 "\n";

        }else  { # if no average, add a \n for each line
            print SF "\n";
            print SF2 "\n";
            print SF3 "\n";
        }
}


##########################################################################

sub load_SpaceDelimFile  {

    # this function loads a table that is space delimited
    # first row is header
    # does not remove any header!

    my $file = shift;
    my @data;

    my $old_separator = $/;
    undef $/;
    open (FILE, $file)  || die "Cannot open $file\n";
    my $file_string = <FILE>;
    close FILE;
    $/ = $old_separator;

    my @tmp = split /\n/, $file_string;
    for (my $r=0; $r<scalar(@tmp);  $r++)  {
	my @row = split /\s/, $tmp[$r];

	for (my $c=0; $c<scalar(@row); $c++)  {
	    $data[$r][$c] = $row[$c];
	}
    }
    my $size = scalar(@data);
    print "file has $size rows\n";

    return @data;
}

##########################################################################
