#!/usr/bin/env perl
#;-*- Perl -*-2

# Set up a dimer saddle point search from a NEB run.
# I am getting the configurations from the CONTCARs.  
# If there is an input argument it will be image number to be used.

use FindBin qw($Bin);
use lib "$Bin";

if($ARGV[0]) {
    $i = $ARGV[0];
    $f = 0;
} else {
    # Select the image from the exts.dat file
    open EXT , "exts.dat"
        or die " No exts.dat in this directory; please run nebresults.pl \n";
    @ext = <EXT>;
    close EXT;
    chomp(@ext);
    $n = @ext;
    $max = -1.0;
    for($i=0; $i<$n; $i++) {
        $t = $ext[$i];
        $t =~ s/^\s+//g;
        @t = split /\s+/,$t;
        if($t[8] > $max) {
            $max = $t[8];
            $im = $t[5];
        }
    }
    $i = int($im);
    $f = $im - $i;
    if($im<=0) {
        $i=0;
        $f=0;
    }
}
#print "\n FORMING DIMER BETWEEN IMAGES ",$i," and ",$i+1,"\n\n";

die "Number of images must be less than 99 for now ! \n" if $i >= 99;

$d1 = sprintf "%02d",$i;
$d2 = sprintf "%02d",$i+1;
open (DIST, "$Bin/dist.pl $d1/POSCAR $d2/POSCAR |");
$d12 = <DIST>;
close DIST;
$d12>0 || die "Zero distance between images ",$i1," and ",$i2,"\n";
$NdR = 5e-3;
$f1 = $f - $NdR/$d12;
$f2 = $f + $NdR/$d12;

# Generate the dimer POSCARs.

system "$Bin/posinterp.pl $d1/CONTCAR $d2/CONTCAR $f1 > /dev/null";
rename "POSCAR.out" , "POSCAR1";
system "$Bin/posinterp.pl $d1/CONTCAR $d2/CONTCAR $f2 > /dev/null";
rename "POSCAR.out" , "POSCAR2";
system "$Bin/posinterp.pl $d1/CONTCAR $d2/CONTCAR $f  > /dev/null";
rename "POSCAR.out" , "POSCAR";


$dir = "dim" ;
mkdir $dir ; # mkdir "$dir/01" ; mkdir "$dir/02";
#system "cp INCAR KPOINTS $dir";
#if(-e "POTCAR") {
#    system "cp POTCAR $dir";
#} elsif(-e "../POTCAR") {
#    system "cp ../POTCAR $dir";
#}
system "$Bin/modemake.pl POSCAR1 POSCAR2 ; mv MODECAR $dir ; mv POSCAR $dir";
#  system "cp POSCAR1 $dir/01/POSCAR" ;
#  system "cp POSCAR2 $dir/02/POSCAR" ;
#print " FOR DIMER, REMEMBER TO SET IN THE INCAR: \n\n";
#print "ICHAIN = 2 \n";
#print "IOPT = 2 \n";
#print "IBRION = 3 \n";
#print "POTIM  = 0.0 \n";
#print "EDIFF  = 1E-7 \n \n";
#print "DdR       (DEFAULT = 5E-3) \n";
#print "DRotMax   (DEFAULT = 4) \n";
#print "DFNMax    (DEFAULT = 1.0) \n";
#print "DFNMin    (DEFAULT = 0.01) \n\n";

unlink glob "POSCAR1";
unlink glob "POSCAR2";

