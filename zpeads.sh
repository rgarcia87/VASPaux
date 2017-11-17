#!/bin/bash
# Zero-point vibrational energies for adsorbed molecules
# Rodrigo García-Muelas
# Dec 2, 2016. Tarragona.
# 
# This program computes the zpe from the OUTCAR file

# Total number of frequencies calculated
free1=`grep "f  =" OUTCAR | wc -l`

# Extracting frequencies to be added
touch frq.dat && rm -f frq.dat
grep "f  =" OUTCAR | head -n $free1 | awk '{print $10}' >> frq.dat

# Adding up energies
d=0 ; sum=0 ; until [ $d -eq $free1 ] ; do 
d=$(($d+1))
add=`sed "${d}q;d" frq.dat` 
sum=`echo "scale=16; 0.0005*$add+$sum " | bc -l`
#echo add $add sum $sum
done

rm -f sumfrq.txt 2>/dev/null
echo $sum >> sumfrq.txt
