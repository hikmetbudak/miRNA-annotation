#!perl -w
# SUmirFind.pl - a script that uses ncbiBLAST to search for potential homologs of known miRNAs
# Written by Stuart J. Lucas, 2012
# Version History:
# v1.1, 23.01.2013: Added option to specify number of mismatches, more feedback during running, and set to delete temporary files.

print "\n\n+++ SUmirFind.pl v1.1 +++\n\n";

my ( $mirnaquery, $blastdatabase ) = @ARGV or die "Please specify a fasta file containing the miRNA sequences you wish to search with, and the BLAST database you wish to search, with full path if it is not in the current directory";

# Get number of mismatches from user

print "Please specify the maximum number of base mismatches you wish to allow between your query miRNAs and test sequences.\n";
print "Maximum permitted mismatches (1-3 recommended): ";
my $mislimit = <STDIN>;
chomp $mislimit;
$mislimit =~ s/\r//;
if ( $mislimit =~ m/\D+/ )
{
die "Your entry was not a number.  Please run again and enter digits only.\n";
}
else
{
print "Searching for miRNA sequences with a maximum of $mislimit mismatches.\n";
}

# Convert fasta file into a table of miRNAs

my $mirnacount = 0;
open ( MIRNAS, $mirnaquery ) or die "Could not open $mirnaquery.  Is it in the right folder?";
open ( MIRNATABLE, ">".$mirnaquery.".tbl" ) or die "Could not open an output file!";

print MIRNATABLE "# miRNA ID\tsequence";

while($line = <MIRNAS>)
{
	chomp $line;	
	$line =~ s/\r//;
	if ($line =~ /^>(\S+)\s*/)
	{
		print MIRNATABLE "\n$1";
		print MIRNATABLE "\t";
		$mirnacount++;
	}
	else
	{
		print MIRNATABLE "$line";
	}
}
print MIRNATABLE "\n";

close MIRNAS;
close MIRNATABLE;
print "$mirnacount miRNA sequences detected in query file.\n";

# Populate hash table of miRNAs from newly generated table

my (%QuerymiRNAs);
open ( MIRNADATA, $mirnaquery.".tbl" ) or die "The miRNA table was not generated";

while ($line = <MIRNADATA>)
{
	chomp $line;
	if ( $line =~ /^#/ )
	{
		next;
	}		
	elsif ( $line =~ /(\S+)\t(\S+)/ )
	{
		my $idkey = $1;
		my $seqval = $2;		
		$QuerymiRNAs{"$idkey"} = "$seqval";
	} 
}

close MIRNADATA;
my $mirnacount2 = scalar(keys(%QuerymiRNAs));

print "miRNA data analysed....  Now running BLAST for $mirnacount2 miRNAs.  This could take some time.\n";

# Run BLAST for the specified miRNAs

system ("blastn -task blastn-short -query $mirnaquery -db $blastdatabase -ungapped -penalty -1 -reward 1 -outfmt 6 -out $mirnaquery.allhits" );

# Filter blast hit table (in output format 6) to remove alignments with >specified mismatches

print "BLAST complete, now filtering.\n";

open ( BLASTHITS, $mirnaquery.".allhits" ) or die "Couldn't find the BLAST output!";
open ( FILTEREDHITS, ">".$mirnaquery.".results.tbl" ) or die "Results file failure";

print FILTEREDHITS "# Query ID\tSubject ID\t%\tlength\tmism\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n";

my $badcount = 0;
my $goodcount = 0;
while ( $line = <BLASTHITS> )
{
	chomp $line;	
	my ( $qid, $sid, $percent, $allength, $mismatch, $gaps, $qstart, $qend, $sstart, $send, $evalue, $bitscore ) = split /\t/, $line;
	my $qlength = length $QuerymiRNAs{ $qid };
	my $difference = $qlength - $allength;
	if ( $mismatch + $difference > $mislimit )
	{
		$badcount++;		
		next;
	}
	else 
	{
		print FILTEREDHITS "$line\n";
		$goodcount++;
	}	
}
close BLASTHITS;
close FILTEREDHITS;
unlink ($mirnaquery.".tbl");
unlink ($mirnaquery.".allhits");

print "Filtering complete.  $badcount hits were rejected due to being too short, or having too many mismatches.\n";
print "$goodcount hits were recorded in the file ".$mirnaquery.".results.tbl\n";

