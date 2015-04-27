#!/usr/bin/perl -w
# merge the info for each CpG sites
# watson and crick Cs for the same CpG dyad are merged based on the CpG location list file
# analyze only one chr at the same time; thus, the cpg.pos and sam.txt file need to be for single chrom.
# Ming-an Sun, May 21, 2013
# Ming-an Sun, Sep 04, 2013
use warnings;
use strict;

=head1 Description

This script is used to parse Bismark result to get the methylation pattern for each CpG dyads.
Both the Bismark result and CpG coordinate need to be provided as input. The information for
the two paired Cs in the same CpG dyad are merged as one based on the CpG coordnate file.

=head1 Usage

perl merge_CpG_site.pl <cpg coordnate file> <bismark output> > output

=head1 Example

perl merge_CpG_site.pl hg19.CPG CpG_context.hg19.fetal.sam.txt > CPG_context.hg19.fetal.sam.txt.site

=head1 Version

Ver 1.1, Ming-an Sun, May 18, 2014 (for multiple chr)

Ver 1.0, Ming-an Sun, Sep 04, 2013 (for individual chr)

=cut

die `pod2text $0` unless @ARGV==2;

my ($cpgPosFile, $samTxtFile) = @ARGV;

# read in CpG position (in + strand)
print STDERR "Parsing $cpgPosFile ...\n";
my %positions;
my $con = 0;
open(POS, $cpgPosFile=~/\.gz$/ ? "zcat -c $cpgPosFile |" : $cpgPosFile)||die"Cannot open $cpgPosFile\n";
while(<POS>){
    if(/^(chr\S+)\s+(\d+)/){
        $con++;
        print "$con\n" if $con%1000000==0;
        $positions{$1}->{$2} = '';
    }
}
close POS;

# parse sam.txt file, and store the information for each CpG dyad
print STDERR "Parsing $samTxtFile ...\n";
my %info;
$con = 0;
open(TXT, $samTxtFile=~/\.gz$/ ? "zcat -c $samTxtFile |" : $samTxtFile)||die"Cannot open $samTxtFile\n";
while(<TXT>){
	next if /^Bismark/;
    $con++;
    print STDERR "$con\n" if $con%1000000==0;
	my ($readName, $pattern, $chr, $position) = split;
	unless(defined $positions{$chr}->{$position} || defined $positions{$chr}->{$position-1}){
        print STDERR "CpG position undefined: $_";
        next;
    }
	$position-- if !defined($positions{$chr}->{$position}) && defined($positions{$chr}->{$position-1}); # adjust the position for crick Cs
	# store the names of reads with mC in [0], and C in [1]
	if(!defined $info{$chr}->{$position}){ # initialize
		if($pattern eq "+"){
			$info{$chr}->{$position}->[0] = "$readName;";
			$info{$chr}->{$position}->[1] = "";
		}else{
			$info{$chr}->{$position}->[0] = "";
			$info{$chr}->{$position}->[1] = "$readName;";
		}
	}else{ # append new info
		if($pattern eq "+"){
			$info{$chr}->{$position}->[0] .= "$readName;";
		}else{
			$info{$chr}->{$position}->[1] .= "$readName;";
		}
	}
}
close TXT;

undef %positions;

# output information for each CpG dyad
print STDERR "Output result ...\n";
foreach my $chr (sort keys %info){
	print STDERR "$chr\n";
	foreach my $pos (sort {$a <=> $b} keys %{$info{$chr}}){
		my $posStr = $info{$chr}->{$pos}->[0];
		my $negStr = $info{$chr}->{$pos}->[1];
        #my $posNum = &getNum($posStr);
        #my $negNum = &getNum($negStr);
		$posStr = 'NA;' if $posStr eq '';
		$negStr = 'NA;' if $negStr eq '';
        print "$chr\t$pos\t$posStr\t$negStr\n";
        #my $allNum = $posNum + $negNum;
        #my $mlev = $allNum>0 ? sprintf("%.5f", $posNum/$allNum) : 'NA';
        #print "$chr\t$pos\t$allNum\t$posNum\t$mlev\t$posStr\t$negStr\n";
	}
    delete $info{$chr};
}

print STDERR "Job finished!\n";

######################################################
################ subroutines #########################

# get number based on the occurrent of ";"
sub getNum{
	my $str = shift;
	my $num = 0;
	while($str =~ /\;/g){
		$num++;
	}
	return $num;
}

__END__
==> chr22.CPG <==
16050097	+
16050098	-
16050114	+
16050115	-
16050174	+
16050175	-
16050206	+
16050207	-
16050213	+
16050214	-

Bismark methylation extractor version v0.7.4
SRR478983.5	+	chr22	33431552	Z
SRR478983.5	+	chr22	33431588	Z
SRR478983.5	+	chr22	33431626	Z
SRR478983.74	+	chr22	45815199	Z
SRR478983.98	+	chr22	20826842	Z
SRR478983.98	+	chr22	20826784	Z
SRR478983.98	+	chr22	20826751	Z
SRR478983.145	-	chr22	45713129	z
SRR478983.204	+	chr22	46758048	Z
