#!/bin/bash
# Create a POTCAR file for a VASP calculation by concatenation of POTCAR files from a repository

# Define pseudopotentials repository:
folder="/home/userxxx/path/to/pp" 

# Check if older version of POTCAR is present
if [ -f POTCAR ] ; then
 mv -f POTCAR old-POTCAR
 echo " ** Warning: old POTCAR file found and renamed to 'old-POTCAR'."
fi 

# Main loop - concatenate the appropriate POTCARs (or archives)
for i in $*
do
 if test -f $repo/$i/POTCAR ; then
  cat $folder/$i/POTCAR >> POTCAR
 elif test -f $repo/$i/POTCAR.Z ; then
  zcat $folder/$i/POTCAR >> POTCAR
 elif test -f $repo/$i/POTCAR.gz ; then
  gunzip -c $folder/$i/POTCAR.gz >> POTCAR
 else
  echo "No suitable POTCAR for element '$i' found!"
 fi
done
