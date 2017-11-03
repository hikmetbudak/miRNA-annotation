#!/usr/bin/perl -w

# SUmirFold.pl - a script that, using a BLAST results table, retrieves sequences from the BLAST database and obtains their predicted secondary structure using UNAfold, after which viable hairpins are detected and retrieved.
# Version 1.0.17, 25/11/14  - Added lines to deal with predictions where putative miRNA starts before the start of the target sequence being searched.  
# 								Such entries are now reported in the results table with initial 'n's for the bases missing from the sequence.
# Version 1.0.16, 26/05/14  - Added lines to avoid crashing if putative hairpin contains only As & Ts
# Version 1.0.15, 30/05/13  - corrected bugs that reported mature miRNA position in hairpin as 1 fewer than the reality, reported incorrect miRNA* sequence
# Version 1.0.14, 15/02/13	- add in code to report miRNA* position and sequence
# Version 1.0.13, 28/01/13  - corrected some missing newlines, changed the filesystem calls to Perl commands to increase portability, delete temporary results files
#                             implemented option to limit the number of results for a single miRNA, and initial check for UNAFold installation
# Version 1.0.12, 17/09/12  - corrected commands for reading in fasta files to deal with Windows carriage returns
# Version 1.0.11, 30/07/12  - retains all of subject ID

my $goodcount = 0;
my $badcount = 0;
my $suspectcount = 0;
my $loopno = 0;

if ( $#ARGV != 2 )
{
	print "Correct usage is:\n";
	print "perl SUmirFold.pl [miRNA query file] [HitTable] [BLASTdb]\n";
}

my ( $mirnaquery, $infile, $blastdb ) = @ARGV or die "Please specify names of a file containing known miRNAs (in fasta format), BLAST hit table (generated using '-outfmt 6' when running BLAST) and blast database from which the hit table was generated.\n";

# Check for presence of UNAFold and maximum number of folds

print "---  SUmirFold.pl v1.0.13  ---\n";
system ("UNAFold.pl -V > unachk.tmp");
open ( UNAREPORT, "unachk.tmp") or die "UNAFold is not available.  Please check that it is installed and in your PATH.\n";
while ($firstline = <UNAREPORT>)
{
	chomp $firstline;
	$firstline = substr ($firstline, 21);
	print "UNAFold version $firstline detected\n";
	last;
}
unlink ("unachk.tmp");
print "miRNAs that produced a large number of hits with SUmirFind\n may be matching repeat elements instead of real miRNAs.\nTo save computing time, you can limit the number of hits \nto be folded for any miRNA.\n";
print "Please enter the maximum no. of hits to be folded (press Return for no limit): ";
my $maxhits = <STDIN>;
chomp $maxhits;
if ( $maxhits )
{
	if ( $maxhits =~ m/\D+/ ) { die "That wasn't a number.  Try again!\n" }
	print "Proceeding to fold a maximum of $maxhits hits per miRNA.\n";
}
else
{
	print "Folding all hits.\n";
}

# Convert fasta file into a table of miRNAs

open ( MIRNAS, $mirnaquery ) or die "Could not open $mirnaquery.  Is it in the right folder?\n";
open ( MIRNATABLE, ">".$mirnaquery.".tbl" ) or die "Could not open an output file!\n";

print MIRNATABLE "# miRNA ID\tsequence";

while ($line = <MIRNAS>)
{
	chomp $line;
	$line =~ s/\r//;
	if ($line =~ /^>(\S+)\s*/)
	{
		print MIRNATABLE "\n$1";
		print MIRNATABLE "\t";
	}
	else
	{
		print MIRNATABLE "$line";
	}
}
print MIRNATABLE "\n";

close MIRNAS;
close MIRNATABLE;

# Populate hash table of miRNAs from newly generated table

my (%QuerymiRNAs);
open ( MIRNADATA, $mirnaquery.".tbl" ) or die "The miRNA table was not generated\n";

while ($line = <MIRNADATA>)
{
	chomp $line;
	if ( $line =~ /^#/ )
	{
		next;
	}		
	elsif ( $line =~ /(\S*)\t(\S*)/ )
	{
		my $idkey = $1;
		my $seqval = $2;		
		$QuerymiRNAs{"$idkey"} = "$seqval";
	} 
}
close MIRNADATA;
unlink ( "$mirnaquery.tbl" );
print "miRNA data analysed.  Preparing to fold RNA sequences.\n";

# Reading in data from BLAST hit table and carrying out initial folds:

unless ( -e $infile && -f $infile && -r $infile )
{
	die "$infile cannot be accessed.  Does it exist?\n";
}

open ( IN, $infile ) or die "I don't have permission to open $infile!\n";
open ( OUT, ">".$infile.".seqtable" ) or die "Can't open an output file!\n";
open ( OUTTABLE, ">".$infile.".tbl") or die "Can't open an output file!\n";
open ( DISCARD, ">".$infile.".rejects" ) or die "Can't open an output file!\n";
open ( SUSPECT, ">".$infile.".suspect.tbl" ) or die "Can't open an output file!\n";
open ( SUSPECTOUT, ">".$infile.".suspect.seqtable" ) or die "Can't open an output file!\n";
open ( LOGFILE, ">".$infile.".log" ) or die "Can't open an output file!\n";
print DISCARD "Matches rejected:\n";
print DISCARD "Seq ID\tmiRNA\tStart\tEnd\tReason for rejection\n";
print OUTTABLE "#Structure       \tMatched    \tConserved\tMatch\tMature miRNA\tComplement\n";
print OUTTABLE "#Filename        \tSequence ID\tmiRNA    \tlength\tstart\tend\tstart\tend\n";
print SUSPECT "#Structure       \tMatched    \tConserved\tMatch\tMature miRNA\tComplement\n";
print SUSPECT "#Filename        \tSequence ID\tmiRNA    \tlength\tstart\tend\tstart\tend\n";

if ( -e "$infile.initialfolds" )
{
	my $nowdate = localtime;
	my $prefix = substr ($nowdate, 4, 12);
	$prefix =~ s/\s|://g;
	my $oldinit = $infile.".initialfolds";
	my $newinit = $prefix."_$oldinit";
	my $oldhairpins = $infile.".hairpins";
	my $newhairpins = $prefix."_$oldhairpins";
	print "A folder called $infile.initialfolds already exists.\n Renaming old results folders with prefix $prefix\n";
	print LOGFILE "A folder called $infile.initialfolds already existed.\n Renamed old results folders with prefix $prefix\n";
	rename ( $oldinit, $newinit );
	rename ( $oldhairpins, $newhairpins );
}
unless ( mkdir "$infile.initialfolds" )  {die "Unable to create results directory\n"};
unless ( mkdir "$infile.hairpins" )  {die "Unable to create results directory\n"};

my $samecount = 0;
my $samename = "";
while ($line = <IN>)
{
	chomp $line;
	$line =~ s/\r//;
	next if $line =~ /^#/;
	next if $line eq "";
	print ".\n";

# Get the values from the table, check whether the maximum folds criterion has been exceeded;
# if not, retrieve the sequence from the BLAST database and write it to a fasta file

	my ( $qid, $sid, $percent, $allength, $mismatch, $gaps, $qstart, $qend, $sstart, $send, $evalue, $bitscore ) = split /\t/, $line;
	$loopno++;
	my $uniqueid = $loopno."_".$sid;
	if ( $maxhits )
	{
		if ( $qid ne $samename )
		{
			$samename = $qid;
			$samecount = 0;
		}
		$samecount++;
		if ( $samecount > $maxhits )
		{
			print LOGFILE "Maximum hits for $qid exceeded, ignoring hit $uniqueid\n";
			next;
		}
	}
	my $qlength = length $QuerymiRNAs{ $qid };	
	system ("blastdbcmd -db $blastdb -entry $sid -outfmt %f -out $uniqueid.fsa" );
	print LOGFILE "Testing hit for $qid on sequence $uniqueid\n";

# Reverse complement the sequence if it is on the negative strand; otherwise just convert Ts to Us.  Either way, clip the sequence if it is too long.
# Also, truncate the header line if it is longer than 130 characters (otherwise UNAFold doesn't generate structure files)

	if ( $sstart > $send )
	{
		print LOGFILE "This hit is on the negative strand.  Initiating reverse complement.\n";		
		my $tempseq = "";
		my $defline = "";
		my $clipstart = 0;		
		open ( NEGSTRAND, $uniqueid.".fsa" );
		while ( $line = <NEGSTRAND> )
		{
			chomp $line;
			if ( $line =~ m/^>/)
			{ 
				$defline = $line;
				if ( length($defline) > 130 )
				{
					$defline = substr ($defline, 0, 130);
				}
			}
			else 
			{
				$tempseq = $tempseq.$line;
			}
		}
		close NEGSTRAND;
		$tempseq = scalar reverse ("$tempseq");
		$tempseq =~ tr/[UT]/a/;
		$tempseq =~ tr/A/u/;
		$tempseq =~ tr/C/g/;
		$tempseq =~ tr/G/c/;
		$tempseq =~ tr/acgu/ACGU/;
		$seqlength = length ($tempseq);
		$sstart = 1 + $seqlength - $sstart;
		$send = 1 + $seqlength - $send;
		if ($seqlength > 750)
		{
			print LOGFILE "Long sequence, clipping 300bp either side of miRNA\n";
			$clipstart = $sstart - 300;
			if ($clipstart < 0)
			{
				$clipstart = 0;
			}
			$tempseq = substr $tempseq, $clipstart, 725;
			$sstart = $sstart - $clipstart;
			$send = $send -$clipstart;
		}
		open  ( POSSTRAND, ">".$uniqueid.".fsa" );		
		print POSSTRAND "$defline\n$tempseq\n";
		close POSSTRAND;			
		print LOGFILE "Reverse complement stats (id, length, hit start and end): $sid, $seqlength, $sstart, $send\n";
	}
	else
	{
		my $tempseq = "";
		my $defline = "";
		open ( DNASTRAND, $uniqueid.".fsa" );
		while ( $line = <DNASTRAND> )
		{
			chomp $line;
			if ( $line =~ m/^>/)
			{ 
				$defline = $line;
				if ( length($defline) > 130 )
				{
					$defline = substr ($defline, 0, 130);
				}
			}
			else 
			{
				$tempseq = $tempseq.$line;
			}
		}
		close DNASTRAND;
		$tempseq =~ tr/T/U/;
		if (length ($tempseq) > 750)
		{
			print LOGFILE "Long sequence, clipping 300bp either side of miRNA\n";
			$clipstart = $sstart - 300;
			if ($clipstart < 0)
			{
				$clipstart = 0;
			}
			$tempseq = substr $tempseq, $clipstart, 725;
			$sstart = $sstart - $clipstart;
			$send = $send - $clipstart;
		}
		open  ( RNASTRAND, ">".$uniqueid.".fsa" );		
		print RNASTRAND "$defline\n$tempseq\n";
		close RNASTRAND;
	}

# Adjust ends of putative mature miRNA if it is shorter than query miRNA, and note 'extended' miRNA duplex end	
	$sstart = $sstart - $qstart + 1;
	$send = $send + $qlength - $qend;
	my $exstart = $sstart - 2; 

# Run UNAfold on the fasta file
	print LOGFILE "Running UNAfold on $sid\n";
	system ("UNAFold.pl -X 1 --ann ss-count $uniqueid.fsa");

# Check that UNAfold worked
	unless (-e $uniqueid.".fsa.plot" && -e $uniqueid.".fsa.ct")
	{
		print "UNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
		print LOGFILE "UNAfold gave no result.  Probably a low complexity sequence. Files deleted.\n";
		print DISCARD "$sid\t$qid\t$sstart\t$send\tNo result from UNAFold\n";
		unlink ( <$uniqueid.*> );
		$badcount++;
		next;
	}
	
# Examine the Ct file for unpaired bases in the miRNA sequence
	open ( FOLD, $uniqueid.".fsa.ct" ) or die "where is the ct file?";
	my $linecount = 0;
	my $ssflag = 0;
	my $revstart = 0;
	my $revstartbackup = 0;
	my $revend = 0;
	my $revendbackup = 0;
	my $revlength = 0;
	while ($line = <FOLD>)
	{
		chomp $line;		
		next if $line =~ /dG/;
		$linecount++;
		next if $linecount < ($exstart);
		next if $linecount > ($send);
		if ( $line =~ /\d*\t\w\t\d*\t\d*\t(\d*)\t\d*\t\d*\t\d*/ )
		{
			$ssflag++ unless $1 > 0; 
			if ( $linecount == $exstart )
			{
				$revend = $1;
			}
			if ( $linecount == $sstart )
			{
				$revendbackup = $1;
			}
			if ( $linecount == $send - 2)
			{
				$revstart = $1;
			}
			if ( $linecount == $send )
			{
				$revstartbackup = $1;
			}
		}
	}
	close FOLD;
	if ( $revend == 0 )
		{
			$revend = $revendbackup - 2;
		}
	if ( $revstart == 0 )
		{
			$revstart = $revstartbackup - 2;
		}
	$revlength = $revend - $revstart + 1;
	if ( $ssflag > 4 )
	{
		print LOGFILE "Too many unpaired bases in the miRNA region of $sid. Files deleted.\n";
		print DISCARD "$sid\t$qid\t$sstart\t$send\t$ssflag unpaired bases in miRNA\n";
		unlink ( <$uniqueid.*> );
		$badcount++;
	}
	elsif ( $revend == -2 or $revstart == -2)
	{
		print LOGFILE "One end of the miRNA region was not base-paired, so could not locate the end of the miRNA complementary sequence.  This hit will be placed in the 'suspect' table.\n";
		if ( $revend == -2 )
		{
			$revend = $revstart + $qlength - 1;
		}
		elsif ( $revstart == -2 )
		{
			$revstart = $revend - $qlength + 1;
		}
		print SUSPECT "$uniqueid\t$sid\t$qid\t$qlength\t$sstart\t$send\t$revstart\t$revend\n";
		open ( FASTA, $uniqueid.".fsa");
		while ($line = <FASTA>)
		{	
			chomp $line;		
			if ( $line =~ /^>/ )
			{
				print SUSPECTOUT "$uniqueid\t";
			}
			else
			{
				print SUSPECTOUT "$line";
			}
		}
		print SUSPECTOUT "\n";
		close FASTA;
		unlink ( <$uniqueid.fsa.*> );
		unlink ( "$uniqueid.fsa_1.ss" );
		unlink ( "$uniqueid.fsa_1.pdf" );
		rename ( "$uniqueid.fsa", "$infile.initialfolds/$uniqueid.fsa" );
		$suspectcount++;
	}	
	elsif ( $revlength - $qlength > 3 )
	{
		print LOGFILE "The miRNA complementary sequence of $sid contains breaks or a large loop. Files deleted.\n";
		print DISCARD "$sid\t$qid\t$sstart\t$send\tmiRNA complementary sequence is broken \n";
		unlink ( <$uniqueid.*> );
		$badcount++;
	}
	elsif ( $ssflag == 0 )
	{
		print LOGFILE "The putative miRNA sequence of $sid is perfectly base-paired, so it is more likely to be an inverted repeat or siRNA.  This hit will be placed in the 'suspect' table.\n";
		print SUSPECT "$uniqueid\t$sid\t$qid\t$qlength\t$sstart\t$send\t$revstart\t$revend\n";
		open ( FASTA, $uniqueid.".fsa");
		while ($line = <FASTA>)
		{	
			chomp $line;		
			if ( $line =~ /^>/ )
			{
				print SUSPECTOUT "$uniqueid\t";
			}
			else
			{
				print SUSPECTOUT "$line";
			}
		}
		print SUSPECTOUT "\n";
		close FASTA;
		unlink ( <$uniqueid.fsa.*> );
		unlink ( "$uniqueid.fsa_1.ss" );
		unlink ( "$uniqueid.fsa_1.pdf" );
		rename ( "$uniqueid.fsa", "$infile.initialfolds/$uniqueid.fsa" );
		$suspectcount++;
	}	
	else
	{
		print LOGFILE "The secondary structure of $uniqueid passes initial analysis, writing to output table and fasta file.\n";		
		print OUTTABLE "$uniqueid\t$sid\t$qid\t$qlength\t$sstart\t$send\t$revstart\t$revend\n";		
		open ( FASTA, $uniqueid.".fsa");
		while ($line = <FASTA>)
		{	
			chomp $line;		
			if ( $line =~ /^>/ )
			{
				print OUT "$uniqueid\t";
			}
			else
			{
				print OUT "$line";
			}
		}
		print OUT "\n";
		close FASTA;
		unlink ( <$uniqueid.fsa.*> );
		unlink ( "$uniqueid.fsa_1.ss" );
		unlink ( "$uniqueid.fsa_1.pdf" );
		rename ( "$uniqueid.fsa", "$infile.initialfolds/$uniqueid.fsa" );
		$goodcount++;
	}
}
close IN;
close OUT;
close OUTTABLE;
close SUSPECT;
close SUSPECTOUT;

# Shunt RNA secondary structures into a folder

my @strucfiles = <*.fsa_1.*>;
foreach $strucfile (@strucfiles)
{
	rename ($strucfile, "$infile.initialfolds/$strucfile");
}

print LOGFILE "From the input table, $goodcount sequence(s) gave folds that could contain a miRNA.\n$badcount sequence(s) were rejected after folding.\n$suspectcount sequence(s) are more likely to be repeats or siRNAs.\n\n";
close LOGFILE;

# Conduct further analysis on hairpin regions of all good hits

open ( GOODTABLE, $infile.".tbl" );
open ( RESULTS, ">".$infile.".hairpins.tbl" ) or die "Can't open an output file!\n";
open ( HAIRPINS, ">".$infile.".hairpins.fsa" ) or die "Can't open an output file!\n";
open ( LOGFILE2, ">".$infile.".hairpins.log" ) or die "Can't open an output file!\n";
print RESULTS "Unique\tNew miRNA\t\t\tConserved miRNA\t\t\tSequence\tMature\tMature\tmiRNA*\tmiRNA*\tmiRNA*\tHairpin\tPre-miRNA stats\n";
print RESULTS "Hit ID\tID\tSequence\tLength\tID\tSequence\tMismatch\tID\tStart\tEnd\tStart\tEnd\tSequence\tlocation\tlength\tMFE\tGC%\tMFEI\tstart\tsequence\n";
my $goodhairpins = 0;

while ( $line = <GOODTABLE> )
{
	chomp $line;
	next if $line =~ /^#/;
	my ( $uniqueid, $sseqid, $mirnaid, $length, $sstart, $send, $revstart, $revend ) = split /\t/, $line;
	my ( $armflag, $seq ) = "";
	my ( $hairpinstart, $hairpinend, $hairpinlength ) = 0;
	print LOGFILE2 "Analysing hit $uniqueid for $mirnaid\n";	
	
# Determine which part of the sequence corresponds to putative hairpin and retrieve it from the sequence table
	if ( $sstart > $revstart )
	{
		$armflag = "3'";
		$hairpinstart = $revstart - 20;
		$hairpinend = $send + 20;
	}
	elsif ( $sstart < $revstart )
	{
		$armflag = "5'";
		$hairpinstart = $sstart - 20;
		$hairpinend = $revend + 20;
	}
	else
	{
		print LOGFILE2 "The miRNA co-ordinates for $uniqueid don't make sense!  Skipping it.\n";
		next;
	}
	
	open ( GOODSEQS, $infile.".seqtable" );
	while ( $line = <GOODSEQS> )
	{
		my ( $fseqid, $fseq) = split /\t/, $line;
		if  ( $fseqid =~ m/$uniqueid/ )
		{
			$seq = $fseq;	
			last;
		}
	}
	close GOODSEQS;
	if ( $hairpinend - $hairpinstart < (2*$length) + 40 )
	{
		print LOGFILE2 "The miRNA region of $sseqid goes round the head of the hairpin; discarded.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmiRNA goes round head of hairpin\n";
		next;
	}
 	if ( $hairpinstart < 1 )
	{
		$hairpinstart = 1;
	}
	if ( $hairpinend > length $seq )
	{
		$hairpinend = length $seq;
	}
	$hairpinlength = $hairpinend - $hairpinstart + 1;
	
	my $hairpinseq = substr $seq, ($hairpinstart-1), $hairpinlength;
	my $matseq = substr $seq, ($sstart-1), $length;	
	my $matstart = $sstart - $hairpinstart+1;
	my $matend = $send - $hairpinstart+1;
	print LOGFILE2 "Hairpin stats (start, end, length):\n";
	print LOGFILE2 "$hairpinstart\t$hairpinend\t$hairpinlength\n";
	print LOGFILE2 "Running UNAFold.pl on the hairpin sequence... ";
	open ( PFASTA, ">".$uniqueid.".hairpin.fsa" ) or die "Can't open an output file!\n";
	print PFASTA ">$sseqid\t$mirnaid\n";
	print PFASTA "$hairpinseq\n";
	close PFASTA;

# Run UNAFold on the hairpin, get useful information from the .ct file

	system ("UNAFold.pl -X 1 --ann ss-count $uniqueid.hairpin.fsa");
	open ( HAIRPINFOLD, $uniqueid.".hairpin.fsa.ct" ) or die "where is the ct file?";
	my $GCcount = 0;
	my ( $starstart, $starend ) = 0;
	my $starseq = "Not defined";
	while ( $line = <HAIRPINFOLD> )
	{
		chomp $line;
		if ( $line =~ m/dG =\s(\S+)/ )
		{		
			$mfe = $1;
		}
		elsif ( $line =~ /(\d*)\t(\w)\t\d*\t\d*\t(\d*)\t\d*\t\d*\t\d*/ )
		{
			if ($1 == $matstart)
			{
				$starend = $3 + 2;
			}
			if ($1 == $matend)
			{
				$starstart = $3 + 2;
			}
			$GCcount++ if $2 =~ /(G|C)/;
		}
	}
	close HAIRPINFOLD;
	if ($starend != 2 && $starstart != 2)
	{
		$starseq = substr $hairpinseq, ($starstart-1), ($starend-$starstart+1);
	}
	
# Process data and print results files

	my $mirnafam = substr ($mirnaid, 3);
	my $conseq = $QuerymiRNAs{"$mirnaid"};
	my $GCcomp = 100*$GCcount/$hairpinlength;
	my $amfe = (0-100*$mfe/$hairpinlength);
	if ( $GCcomp == 0 )
	{
		print LOGFILE2 "GC content appears to be 0.  That's weird.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
		unlink ( <$uniqueid.*> );
		next;
	}
	my $mfei = $amfe/$GCcomp;
	my $mism = 0;
	for ($loopcount = 0; $loopcount < $length; $loopcount++)
	{
		$mism++ if substr ($conseq, $loopcount, 1) ne substr ($matseq, $loopcount, 1);
	}
	if ( $GCcomp < 24 or $GCcomp > 71 )
	{
		print LOGFILE2 "The hairpin for $uniqueid has too low or high GC content.  Deleting files.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
		unlink ( <$uniqueid.*> );
	}
	elsif ( $mfei < 0.67 )
	{
		print LOGFILE2 "The hairpin for $uniqueid has MFEI < 0.67.  Deleting files.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tMFEI outside acceptable range\n";
		unlink ( <$uniqueid.*> );
	}
	else
	{
		print LOGFILE2 "The hairpin passes analysis, printing to results table\n";
		print RESULTS "$uniqueid\tcandidate$mirnafam\t$matseq\t$length\t$mirnaid\t$conseq\t$mism\t$sseqid\t$matstart\t$matend\t$starstart\t$starend\t$starseq\t$armflag\t$hairpinlength\t$mfe\t$GCcomp\t$mfei\t$hairpinstart\t$hairpinseq\n";
		print HAIRPINS ">$sseqid\t$mirnaid\n$hairpinseq\n";
		unlink ( <$uniqueid.hairpin.fsa.*> );
		unlink ( "$uniqueid.hairpin.fsa_1.ss" );
		unlink ( "$uniqueid.hairpin.fsa_1.pdf" );
		rename ( "$uniqueid.hairpin.fsa", "$infile.hairpins/$uniqueid.hairpin.fsa" );
		$goodhairpins++
	}
}
close RESULTS;
close HAIRPINS;
close GOODTABLE;

# Carry out the same hairpin analysis on suspect hits
print LOGFILE2 "Moving on to examine suspect hits from initial folding analysis...\n";
open ( SUSPECTTABLE, $infile.".suspect.tbl" );
open ( SUSPECTRESULTS, ">".$infile.".suspecthairpins.tbl" ) or die "Can't open an output file!\n";
open ( SUSPECTHAIRPINS, ">".$infile.".suspecthairpins.fsa" ) or die "Can't open an output file!\n";
print SUSPECTRESULTS "Unique\tNew miRNA\t\t\tConserved miRNA\t\t\tSequence\tMature\tMature\tmiRNA*\tmiRNA*\tmiRNA*\tHairpin\tPre-miRNA stats\n";
print SUSPECTRESULTS "Hit ID\tID\tSequence\tLength\tID\tSequence\tMismatch\tID\tStart\tEnd\tStart\tEnd\tSequence\tlocation\tlength\tMFE\tGC%\tMFEI\tstart\tsequence\n";
my $suspecthairpins = 0;

while ( $line = <SUSPECTTABLE> )
{
	chomp $line;
	next if $line =~ /^#/;
	my ( $uniqueid, $sseqid, $mirnaid, $length, $sstart, $send, $revstart, $revend ) = split /\t/, $line;
	my ( $armflag, $seq ) = "";
	my $unknowns = "STR"; 
	my ( $hairpinstart, $hairpinend, $hairpinlength, $unknownbases ) = 0;
	print LOGFILE2 "Analysing hit $uniqueid\n";	
	
# Determine which part of the sequence corresponds to putative hairpin and retrieve it from the sequence table
	if ( $sstart > $revstart )
	{
		$armflag = "3'";
		$hairpinstart = $revstart - 20;
		$hairpinend = $send + 20;
	}
	elsif ( $sstart < $revstart )
	{
		$armflag = "5'";
		$hairpinstart = $sstart - 20;
		$hairpinend = $revend + 20;
	}
	else
	{
		print LOGFILE2 "The miRNA co-ordinates for $uniqueid don't make sense!  Skipping it.\n";
		next;
	}
	
	open ( SUSPECTSEQS, $infile.".suspect.seqtable" );
	while ( $line = <SUSPECTSEQS> )
	{
		my ( $fseqid, $fseq) = split /\t/, $line;
		if  ( $fseqid =~ m/$uniqueid/ )
		{
			$seq = $fseq;
			last;	
		}
	}
	close SUSPECTSEQS;
	if ( $hairpinend - $hairpinstart < (2*$length) + 40 )
	{
		print LOGFILE2 "The miRNA region of $sseqid goes round the head of the hairpin; discarded.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tmiRNA goes round head of hairpin\n";
		next;
	}
 	if ( $hairpinstart < 1 )
	{
		$hairpinstart = 1;
	}
	if ( $hairpinend > length $seq )
	{
		$hairpinend = length $seq;
	}
	$hairpinlength = $hairpinend - $hairpinstart + 1;
	if ( $sstart < 1 )
	{
		$length -= (1 - $sstart);
		$unknownbases = 1 - $sstart;
		$unknowns = ("n" x $unknownbases);
		$sstart = 1;
		print LOGFILE2 "Putative mature miRNA starts before beginning of hairpin sequence, adjusted start and length.\n";
	}
	my $hairpinseq = substr $seq, ($hairpinstart-1), $hairpinlength;
	my $matseq = substr $seq, ($sstart-1), $length;
	if ( $unknowns ne "STR" )
	{
		$matseq = $unknowns.$matseq;
		$unknowns = "STR";
	}	
	my $matstart = $sstart - $hairpinstart+1;
	my $matend = $send - $hairpinstart+1; 
	print LOGFILE2 "Hairpin stats (start, end, length):\n";
	print LOGFILE2 "$hairpinstart\t$hairpinend\t$hairpinlength\n";
	print LOGFILE2 "Running UNAFold.pl on the hairpin sequence... ";
	open ( PFASTA, ">".$uniqueid.".hairpin.fsa" ) or die "Can't open an output file!\n";
	print PFASTA ">$sseqid\t$mirnaid\n";
	print PFASTA "$hairpinseq\n";
	close PFASTA;

# Run UNAFold on the hairpin, get useful information from the .ct file

	system ("UNAFold.pl -X 1 --ann ss-count $uniqueid.hairpin.fsa");
	open ( HAIRPINFOLD, $uniqueid.".hairpin.fsa.ct" ) or die "where is the ct file?";
	my $GCcount = 0;
	my ( $starstart, $starend ) = 0;
	my $starseq = "Not defined";
	while ( $line = <HAIRPINFOLD> )
	{
		chomp $line;
		if ( $line =~ m/dG =\s(\S+)/ )
		{		
			$mfe = $1;
		}
		elsif ( $line =~ /(\d*)\t(\w)\t\d*\t\d*\t(\d*)\t\d*\t\d*\t\d*/ )
		{
			if ($1 == $matstart)
			{
				$starend = $3 + 2;
			}
			if ($1 == $matend)
			{
				$starstart = $3 + 2;
			}
			$GCcount++ if $2 =~ /(G|C)/;
		}
	}
	close HAIRPINFOLD;
	if ($starend != 2 && $starstart != 2)
	{
		$starseq = substr $hairpinseq, ($starstart-1), ($starend-$starstart+1);
	}

# Process data and print results files

	my $mirnafam = substr ($mirnaid, 3);
	my $conseq = $QuerymiRNAs{"$mirnaid"};
	my $GCcomp = 100*$GCcount/$hairpinlength;
	if ( $GCcomp == 0 )
	{
		print LOGFILE2 "GC content appears to be 0.  That's weird.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
		unlink ( <$uniqueid.*> );
		next;
	}
	my $amfe = (0-100*$mfe/$hairpinlength);
	my $mfei = $amfe/$GCcomp;
	my $mism = 0;
	for ($loopcount = 0; $loopcount < $length; $loopcount++)
	{
		$mism++ if substr ($conseq, $loopcount, 1) ne substr ($matseq, $loopcount, 1);
	}
	if ( $GCcomp < 24 or $GCcomp > 71 )
	{
		print LOGFILE2 "The hairpin for $uniqueid has too low or high GC content.  Deleting files.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tGC content outside acceptable range\n";
		unlink ( <$uniqueid.*> );
	}
	elsif ( $mfei < 0.67 )
	{
		print LOGFILE2 "The hairpin for $uniqueid has MFEI < 0.67.  Deleting files.\n";
		print DISCARD "$sseqid\t$mirnaid\t$sstart\t$send\tMFEI outside acceptable range\n";
		unlink ( <$uniqueid.*> );
	}
	else
	{
		print LOGFILE2 "The hairpin passes analysis, printing to suspect results table\n";
		print SUSPECTRESULTS "$uniqueid\tcandidate$mirnafam\t$matseq\t$length\t$mirnaid\t$conseq\t$mism\t$sseqid\t$matstart\t$matend\t$starstart\t$starend\t$starseq\t$armflag\t$hairpinlength\t$mfe\t$GCcomp\t$mfei\t$hairpinstart\t$hairpinseq\n";
		print SUSPECTHAIRPINS ">$sseqid\t$mirnaid\n$hairpinseq\n";
		unlink ( <$uniqueid.hairpin.fsa.*> );
		unlink ( "$uniqueid.hairpin.fsa_1.ss" );
		unlink ( "$uniqueid.hairpin.fsa_1.pdf" );
		rename ( "$uniqueid.hairpin.fsa", "$infile.hairpins/$uniqueid.hairpin.fsa" );
		$suspecthairpins++
	}
}


print LOGFILE2 "\nFrom $goodcount viable folds, $goodhairpins passed the hairpin statistics criteria.\n";
print LOGFILE2 "From $suspectcount suspect but possible folds, $suspecthairpins passed the hairpin statistics criteria.\n";
my @strucfiles2 = <*.fsa_1.*>;
foreach $strucfile (@strucfiles2)
{
	rename ($strucfile, "$infile.hairpins/$strucfile");
}
close SUSPECTTABLE;
close SUSPECTHAIRPINS;
close SUSPECTRESULTS;
close DISCARD;
close LOGFILE2;
unlink ( "$infile.seqtable" );
unlink ( "$infile.suspect.seqtable" );
