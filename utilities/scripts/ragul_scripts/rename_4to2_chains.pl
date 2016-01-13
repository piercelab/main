#!/usr/bin/perl

use strict;

my $pdb = $ARGV[0];
if ($pdb eq "") { die("usage: rename_4to2_chains.pl pdb \n"); }


# open the pdb file
open(PDB, $pdb) || die("unable to open file: $pdb\n");
my @pdb_lines = <PDB>;
close(PDB);

use File::Basename;
my $filename = basename($pdb,  ".pdb_aho.pdb");
print "$filename\n";

open(TCRAHOPDB,">$filename\_aho_chains.pdb") or die("Could not open file for writing: $!\n");

my $chain = "A";
foreach my $line (@pdb_lines) {
    if ( (substr($line, 0, 4) eq "ATOM") &&  (substr($line, 21, 1) eq $chain) ) {
	print TCRAHOPDB "$line";
    }
}
print TCRAHOPDB "TER\n";

my $chain = "C";
my $newchainid = "A";
foreach my $line (@pdb_lines) {
    if ( (substr($line, 0, 4) eq "ATOM") &&  (substr($line, 21, 1) eq $chain) ) {
	my $ori_resnum = substr($line, 22, 5);
	my $resnum = $ori_resnum+500;
	my $tmpline = $line;
	my $formated_res_num = sprintf("%4s ", $resnum);#residue insertion code substituted with space
	my $formated_chain_id = sprintf("%1s", $newchainid);#residue insertion code substituted with space
	substr($tmpline, 22, 5) = $formated_res_num;
	substr($tmpline, 21, 1) = $formated_chain_id;
	print TCRAHOPDB "$tmpline";
    }
}
print TCRAHOPDB "TER\n";

my $chain = "D";
foreach my $line (@pdb_lines) {
    if ( (substr($line, 0, 4) eq "ATOM") &&  (substr($line, 21, 1) eq $chain) ) {
	print TCRAHOPDB "$line";
    }
}
print TCRAHOPDB "TER\n";

my $chain = "E";
my $newchainid = "D";
foreach my $line (@pdb_lines) {
    if ( (substr($line, 0, 4) eq "ATOM") &&  (substr($line, 21, 1) eq $chain) ) {
	my $ori_resnum = substr($line, 22, 5);
	my $resnum = $ori_resnum+500;
	my $tmpline = $line;
	my $formated_res_num = sprintf("%4s ", $resnum);#residue insertion code substituted with space
	my $formated_chain_id = sprintf("%1s", $newchainid);#residue insertion code substituted with space
	substr($tmpline, 22, 5) = $formated_res_num;
	substr($tmpline, 21, 1) = $formated_chain_id;
	print TCRAHOPDB "$tmpline";
    }
}
print TCRAHOPDB "TER\n";


close(TCRAHOPDB);
