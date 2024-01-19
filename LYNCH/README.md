# Title: Proteomics analysis of Lynch Syndrome patients microbiome
# Last update: 10th of January, 2024
# PI: Prof. Manuel Ferrer
################################################################

# Index
- 1. [Summary](#1.-Summary)
- 1.1. [Aims of the study](#1.1.-Aims-of-the-study)
- 1.2. [Goals for the analysis](#1.2.-Goals-for-the-analysis)
- 2. [Preprocessing of the raw data](#2.-Preprocessing-of-the-raw-data)

# 1. Summary
Lynch syndrome (LS) is an inherited condition involving a high risk of colorectal cancer (CRC), endometrial cancer and other tumors. The increased cancer risk is associated to the accumulation of genetic mutations due to germline mutations in the mismatch repair genes MLH1, MSH2, MSH6 or PMS2. Recent investigations of 16S rRNA gene sequences in colon biopsies and stool samples revealed microbiome changes associated to specific mismatch repair genes and occurred prior to colorectal neoplasia; however, they were not predictors of baseline or interval adenoma occurrence. Here, we investigated changes of the gut microbial proteomes associated with the mismatch repair gene and high genetic risk of adenomatous polyposis and CRC in stool samples of individuals with LS.

We recruited six study participants (50% female; 50% male) at one clinical site in Spain undergoing mismatch repair gene screening (MLH1, 2; MSH2, 2, MSH6, 2) to confirm LS. We extracted the bacterial proteins from stool samples obtained across 3 (n=4), 5 (n=1) and 6 (n=1) time periods ranging from 84 to 494 days (25 samples in total), and we performed mass spectrometry. During the course of the study, three colon polyps and a CRC were detected in one participant and one colon polyp in another.

## 1.1. Aims of the study
We have three different objectives:
- To detect different clustering among samples regarding the different mutations (MLH1, MSH2 or MSH6).
- To check differences among the time serie proteomes regarding the moment when cancer/polyp is detected.
- Check taxonomic differences among proteomes regarding mutations or the cancer presence.
## 1.2. Goals for the analysis
- Gather all proteomes information in one table
- Add proteins annotation
- Calculate the NSAF and relative parameters and combine with emPAI
- Ordination analysis (PCA or NMDS) regarding proteome profile
- Ordination analysis (PCA or NMDS) regarding proteome taxonomy abundance
- Volcano plots
- Quantile plots

# 2. Preprocessing of the raw data
Each individual sample table was copied into .tsv file for their further merging in a unique table. To this general table we will also add the each protein annotation.
## 2.1. Merging sample tables and adding protein annotation
In order to have only one working table with all the information, we will process all samples tables in .tsv format (e.g. ML002-T1_20191111.tsv) and the annotation table (annotationTable.tsv) by a Perl script:
```Perl
#!/usr/bin/perl
use strict;
use Tie::IxHash;
use Sort::Versions;
# ####################################################################################################

my %hash;
my %hash2;
my %hash3;
my %hash4;
tie %hash, "Tie::IxHash";
tie %hash2, "Tie::IxHash";
tie %hash3, "Tie::IxHash";
tie %hash4, "Tie::IxHash";

# Processing each sample file
opendir(DIR,$ARGV[0]) || suddie("can't find directory $ARGV[0]");
      while(my $f=readdir(DIR)){
            chomp($f);
            if($f=~/(.*)_\d+\.tsv$/){
                  my $sam=$1;
                  $hash2{$sam}="";
                  open(TB,"<$ARGV[0]/$f") || die("can't find $ARGV[0]/$f");
                        print "# reading $ARGV[0]/$f\n";
                        while(my $l=<TB>){
                              chomp($l);
                              my @a=split("\t",$l);
                              $hash3{$a[0]}="$a[2]\t$a[3]";
                              $hash{$a[0]}{$sam}{"PSMs"}=$a[5];
                              $hash{$a[0]}{$sam}{"emPAI"}=$a[7];
                        }
                  close(TB);
            }else{
                  next;
            }
      }
closedir(DIR);

# Gathering each protein annotation
 print "# Processing annotation file";
open(ANN,"<$ARGV[1]") || subdie($ARGV[1]);
      while(my $l=<ANN>){
            chomp($l);
            my @b=split("\t",$l);
            $hash4{$b[0]}="$b[2]\t$b[4]\t$b[7]\t$b[8]\t$b[9]\t$b[10]";
      }
close(ANN);

# Creating the table with all samples information and annotations
print "# Creating table\n";
open(OUT,">$ARGV[2]"); # Creating output file
# Header
print OUT"Lynch_ID\tMW (Da)\tpI\tGene_Name\tGene_length\tPhylum\tGenus\tKO\teggNOG\t".join("\t\t",sort {versioncmp($a,$b)} keys %hash2)."\n";
# Body
foreach my $prot(keys %hash3){
      chomp($prot);
      print OUT"$prot\t$hash3{$prot}";
      print "# Searching annotation for $prot\n";
      print OUT"\t$hash4{$prot}";
      print "# Adding samples data for $prot\n";
      foreach my $sample(sort {versioncmp($a,$b)} keys %hash2){
            if(exists $hash{$prot}{$sample}){
                  print OUT"\t$hash{$prot}{$sample}{PSMs}\t$hash{$prot}{$sample}{emPAI}";
            }else{
                  print OUT"\t0\t0";
            }
      }
      print OUT"\n";
}
close(OUT);
exit;
```
