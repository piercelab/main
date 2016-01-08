#!/usr/bin/perl

# taking the structure as input,
# match its sequence to alpha or beta hmm and get the CDR sequences

use strict;

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

my $pdb = $ARGV[0];
my $alphabeta = $ARGV[1];
my $chain = $ARGV[2];

if ($alphabeta eq "") { die("usage: show_cdr_seqs.pl pdb alpha/beta [chain]\n"); }

# open the pdb file
open(PDB, $pdb) || die("unable to open file: $pdb\n");
my @pdb_lines = <PDB>;
close(PDB);

my @residues = (); # contains the mapping of position number to residue ID
my $seq = "";
my $i = 0;
my %res_counted = (); # avoid double-counting for double occupancy CA atoms
foreach my $line (@pdb_lines)
{
    if ((substr($line, 0, 4) eq "ATOM") && (substr($line, 12, 4) eq " CA "))
    {
	if (($chain ne "") && (substr($line, 21, 1) ne $chain)) { next; }
	my $resid = substr($line, 21, 6);
	if ($res_counted{$resid} ne "") { next; }
	my $res = substr($line, 17, 3);
	$seq .= $aa{$res};
	$residues[++$i] = $resid;
	$res_counted{$resid} = 1;
    }
}

my $aa = "ACDEFGHIKLMNPQRSTVWY";
my $regex_alpha = "^[$aa]+C[$aa]([$aa]{8,12}W)[YF][$aa]{13}([$aa]{6,11})[$aa]{2}[RGVMH][$aa]{4}([$aa]{7,9})[$aa]{9}[DL][$aa]{2,3}Y[$aa][CW][$aa]([$aa]{7,16}[FW])G[$aa]G";

my $regex_beta = "^[$aa]+C[$aa]([$aa]{8,12}W)[Y][$aa]{13}([$aa]{6,11})[$aa]{2}[RGVMH][$aa]{4}([$aa]{7,9})[$aa]{9}[YL][$aa]{2,3}Y[$aa][CW][$aa]([$aa]{7,16}[F])G[$aa]G";

my $regex = $regex_alpha;
if ($alphabeta eq "beta") { $regex = $regex_beta; }

if ($seq =~ /$regex/) 
{ 
    print "$pdb\t$1\t$2\t$3\t$4\n"; 
}
else
{
    print "$pdb\tfailed\n";
}




