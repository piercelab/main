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
my @resmap = (); 
foreach my $line (@pdb_lines)
{
    if ((substr($line, 0, 4) eq "ATOM") && (substr($line, 12, 4) eq " CA "))
    {
	if (($chain ne "") && (substr($line, 21, 1) ne $chain)) { next; }
	my $resid = substr($line, 22, 4);
	if ($res_counted{$resid} ne "") { next; }
	my $res = substr($line, 17, 3);
	$seq .= $aa{$res};
	$residues[$i++] = $resid;
	$res_counted{$resid} = 1;
	my $res_seq_num = substr($line, 22, 4);
        push(@resmap, [$aa{$res}, $resid ]);
    }
}
print "SEQ : $seq\n";
my $aa = "ACDEFGHIKLMNPQRSTVWY";
my $regex_alpha = "^[$aa]+C[$aa]([$aa]{8,12}W)[YF][$aa]{13}([$aa]{6,11})[$aa]{2}[RGVMH][$aa]{4}([$aa]{7,9})[$aa]{9}[DL][$aa]{2,3}Y[$aa][CW][$aa]([$aa]{7,16}[FW])G[$aa]G";

my $regex_beta = "^[$aa]+C[$aa]([$aa]{8,12}W)[Y][$aa]{13}([$aa]{6,11})[$aa]{2}[VILG][$aa]{7}([$aa]{7,9})[$aa]{13}+[YL][$aa]+[CW][$aa]([$aa]{7,16}[F])G[$aa]G";

my $regex;
my $chain_id;
if ($alphabeta eq "alpha") {
    $regex = $regex_alpha;
    $chain_id = 'L';
}
elsif ($alphabeta eq "beta") {
    $regex = $regex_beta;
    $chain_id = 'H';
}

if ($seq =~ /$regex/) 
{ 
    print "$pdb\t$1\t$2\t$3\t$4\n"; 
}
else
{
    print "$pdb\t Pattern match failed\n";
}

my $sequence = chomp($seq);

my $outanarci = "/Users/ragul/Desktop/tmp/num.out";
system('/Users/ragul/anarci-1.0/bin/ANARCI', '-o', $outanarci, '-i', $seq, '-r', 'tr', '-s', 'a');

if ( $? == -1 )
{
  print "ANARCI command failed: $!\n";
}


open(my $fh, '<:encoding(UTF-8)', $outanarci)
  or die "Could not open file '$outanarci' $!";

my $newseq = ""; 
my @numarray = ();
while (my $row = <$fh>) {
    $row =~ /^[#\/]/ and next; #ignore lines starts with # and /
    my @chunks = split ' ', $row;
    #print "$chunks[0] $chunks[1] $chunks[2]\n";
    $newseq .= $chunks[2];
    push (@numarray, \@chunks);
}


if (length($seq) == length($newseq))
{
	print "seq match\n";
}
else
{
	print "seq not match\n";


my $regex_seqmatch = "([$aa+])$newseq([$aa+])";
if ($seq =~ /([$aa]*)$newseq([$aa]*)/) 
{
    my $flen = length($1);
	while ($flen > 0)
    {
	my $chain_id = 'L';
	my $residue_id = substr($1, $flen-1, 1);
	my $residue_num = $numarray[0][1] - 1;
	unshift(@numarray, [$chain_id,$residue_num,$residue_id ]); 
	$flen =  $flen - 1;
    }

    my $llen = length($2);
	while ($llen > 0)
    {
	my $residue_id = substr($2, length($2) - $llen, 1);
	my $residue_num = $numarray[-1][1] + 1;
	push(@numarray, [$chain_id,$residue_num,$residue_id ]); 
	$llen =  $llen - 1;
    }

}

}


open(TCRAHOPDB,">$pdb\_aho.pdb") or die("Could not open file for writing: $!\n");

my $resmap;
my $numarray;
my $arrSize1 = @numarray;
my $arrSize2 = @resmap;
if ($arrSize1 == $arrSize2){
    my $counter = 0;
    while($counter < $arrSize1){
	foreach my $line (@pdb_lines)
	{
	    if (substr($line, 0, 4) eq "ATOM")
	    {
		my $res = substr($line, 17, 3);
		my $seq = $aa{$res};
		my $resid = substr($line, 22, 4);
		my $oldchainid = substr($line,21,1);
		if ( ($seq eq $numarray[$counter][2]) && ($resid == $resmap[$counter][1])) {
		    my $formated_res_num = sprintf("%4s", $numarray[$counter][1]);
		    $line =~ s/$resmap[$counter][1]/$formated_res_num/g;
		    substr($line, 21, 1) = $chain_id;#replace chain id
		    print TCRAHOPDB "$line";
		}
	    }
	    
	}    
	
	$counter++;
	
    }
}
close(TCRAHOFILE);
