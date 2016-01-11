#!/usr/bin/perl

# taking the structure as input,
# use ANARCCI to get AHO numbering and renumber the input structure

use strict;

my $pdb = $ARGV[0];
if ($pdb eq "") { die("usage: renumber_pdb_aho.pl pdb \n"); }


my %aa = ('ALA' => 'A',
	  'CYS' => 'C',
	  'ASP' => 'D',
	  'GLU' => 'E',
	  'PHE' => 'F',
	  'GLY' => 'G',
	  'HIS' => 'H',
	  'ILE' => 'I',
	  'LYS' => 'K',
	  'LEU' => 'L',
	  'MET' => 'M',
	  'ASN' => 'N',
	  'PRO' => 'P',
	  'GLN' => 'Q',
	  'ARG' => 'R',
	  'SER' => 'S',
	  'THR' => 'T',
	  'VAL' => 'V',
	  'TRP' => 'W',
	  'TYR' => 'Y');

my $amino = "ACDEFGHIKLMNPQRSTVWY";


# open the pdb file
open(PDB, $pdb) || die("unable to open file: $pdb\n");
my @pdb_lines = <PDB>;
close(PDB);

use List::MoreUtils qw(uniq);
my @chain_id = ();
foreach my $line (@pdb_lines) {
    if (substr($line, 0, 4) eq "ATOM") {
	push(@chain_id, substr($line, 21, 1));
    }
}

my @unique_chain_id = uniq @chain_id;

open(TCRAHOPDB,">$pdb\_aho.pdb") or die("Could not open file for writing: $!\n");

foreach my $chain (@unique_chain_id) { 
    my $seq = "";
    my $res_num = 0;
    my @resmap = (); 
    foreach my $line (@pdb_lines) {
	if ((substr($line, 0, 4) eq "ATOM") && (substr($line, 12, 4) eq " CA ")) {
	    if (($chain ne "") && (substr($line, 21, 1) ne $chain)) { next; }
	    my $res = substr($line, 17, 3);
	    my $resid = substr($line, 22, 5);#Residue sequence number including any insertion code
	    $seq .= $aa{$res};
	    push(@resmap, [$aa{$res}, $resid ]);
	}
    }
    
    my $outanarci = "anarcci.out";
    system('/Users/ragul/Downloads/Re__anarci/anarci-1.0.BP/bin/ANARCI', '-o', $outanarci, '-i', $seq, '-r', 'tr', '-s', 'a');
    if ( $? == -1 ) {
	print "ANARCI command failed: $!\n";
	exit;	
    }
    
    open(my $fh, '<:encoding(UTF-8)', $outanarci) or die "Could not open file '$outanarci' $!";
    my $newseq = ""; 
    my @numarray = ();
    while (my $row = <$fh>) {
	$row =~ /^[#\/]/ and next; #ignore lines starts with # and /
	my @chunks = split ' ', $row;
	#print "$chunks[0] $chunks[1] $chunks[2]\n";
	next if ($chunks[2] eq "-"); #skip gaps
	$newseq .= $chunks[2];
	push (@numarray, \@chunks);
    }
    
    print "SEQ : $seq\n";
    print "nSEQ : $newseq\n";
    if ($newseq eq "") {
	if ( ($chain eq "D") || ($chain eq "E") ) { 
	    print "Chain $chain : No TCR domain detected\n";
	}
	foreach my $line (@pdb_lines) {
	    if ( (substr($line, 0, 4) eq "ATOM") &&  (substr($line, 21, 1) eq $chain) ) {
		print TCRAHOPDB "$line";
	    }
	}
    } else {
#	print "Chain $chain : TCR domain detected\n";
	unless (length($seq) == length($newseq))
	{
	    my $regex_seqmatch = "([$amino+])$newseq([$amino+])";
	    if ($seq =~ /([$amino]*)$newseq([$amino]*)/) 
	    {
		my $flen = length($1);
		while ($flen > 0)
		{
		    my $residue_id = substr($1, $flen-1, 1);
		    my $residue_num = $numarray[0][1] - 1;
		    unshift(@numarray, [$chain,$residue_num,$residue_id ]); 
		    $flen =  $flen - 1;
		}
		
		my $llen = length($2);
		while ($llen > 0)
		{
		    my $residue_id = substr($2, length($2) - $llen, 1);
		    my $residue_num = $numarray[-1][1] + 1;
		    push(@numarray, [$chain,$residue_num,$residue_id ]); 
		    $llen =  $llen - 1;
		}
	    }
	}
	
	
	
	my $resmap;
	my $numarray;
	my $arrSize1 = @numarray;
	my $arrSize2 = @resmap;
	
	for my $i (0 .. $arrSize1) {
	    foreach my $line (@pdb_lines) {
		if ( (substr($line, 0, 4) eq "ATOM") &&  (substr($line, 21, 1) eq $chain) ) {
		    my $res = substr($line, 17, 3);
		    my $seq = $aa{$res};
		    my $resid = substr($line, 22, 5);
		    if ( ($seq eq $numarray[$i][2]) && ($resid == $resmap[$i][1]) ) {
			my $tmpline = $line;
			my $formated_res_num = sprintf("%4s", $numarray[$i][1]);
			$tmpline =~ s/$resmap[$i][1]/$formated_res_num/g;
			print TCRAHOPDB "$tmpline";
		    }
		}
	    }
	}
    }
}
close(TCRAHOFILE);
