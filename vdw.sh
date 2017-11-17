#!/bin/bash
# Rodrigo García-Muelas
# June 3rd, 2016 
#
# DESCRIPTION
# This program check the list of species contained in a given file (i.e. POSCAR), 
# and generates the corresponding POTCAR, and put the right vdW values on INCAR.
# 
# Acknowledgments: To Martin Gumbau for his help debugging the 4th part. 
# 
# INPUT
# $1 Input file  (i.e. POSCAR)
# $2 (empty)                 
# 
# Indexes
# a number of line
# b number of column
# c number of lattice vector, related to a
# d number of element, related to a
#
# VECTORS
# element 
# numelem

# # # - 1 - R E A D   S P E C I E S   F R O M   P O S C A R 

a=0
totalelem=0
totalatom=0
unset element numelem coord totalelem totalatom 
while read linea ; do   ##### 1 #####
a=$(($a+1)) ; b=0 ; d=$(($a-9)) ; line[$a]=$linea ;
# echo $a; #echo $line;  if [ $a -lt 10 ] ; then  echo "${line[$a]}" ;  fi
for word in ${line[$a]} ; do   ##### 2 #####
b=$(($b+1))
#if [ $a -lt 10 ] ; then echo "previous to case a is $a" ; fi
case $a in
#1) title=${line[$a]} ;;
#2) scaling=$word ;;
#[3-5])  c=$(($a-2)) ; lattice[$c$b]=$word ;; #; echo "coordinates $c $b $word ${lattice[$c$b]}" 
6) totalelem=$(($totalelem+1)) ; element[$b]=$word ;;
7) numelem[$b]=$word ; totalatom=$(($totalatom+$word))  ;;
#8) selective=${line[$a]} ;;
#9) directcartesian=${line[$a]} ;;
#*) coord[$d$b]=$word ;; # echo "${coord[$d$b]}" ;;
esac
#if [ $a -lt 10 ] ; then echo "word value is $word" ; fi
done ##### 2 #####
done < "$1"  ##### 1 #####

# # # - 2 - G E N E R A T E   P O T C A R 

~/bin/potgen ${element[*]}

# # # - 3 - C L E A N   I N C A R 

sed -i '/van der Waals/d' INCAR
sed -i '/LVDW/d'          INCAR
sed -i '/VDW_VERSION/d'   INCAR
sed -i '/Parameters/d'    INCAR
sed -i '/VDW_C6/d'        INCAR
sed -i '/VDW_R0/d'        INCAR

# # # - 4 - G E N E R A T E   V D W   P A R T 

#echo   "element"
#echo   ${element[*]}
function join { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }
list1=`join , ${element[*]} `;
#echo "this is list1: $list1" 

# This option do not recognize arguments in the list with more than one character
#echo "paste   ~/bin/vdw-data/[VDW-TMP.tMp,$EEEE] "
#paste   ~/bin/vdw-data/[VDW-TMP.tMp,$EEEE] 

echo "#!/bin/bash "  >> VDW-TMP.eXe
if [ -z ${element[2]} ] ; then 
   echo "paste -d \" \" ~/bin/vdw-data/VDW-TMP.tMp ~/bin/vdw-data/${element[1]} >> INCAR" >> VDW-TMP.eXe
   else
   echo "paste -d \" \" ~/bin/vdw-data/VDW-TMP.tMp ~/bin/vdw-data/{$list1}      >> INCAR" >> VDW-TMP.eXe
   fi
chmod +x VDW-TMP.eXe ; ./VDW-TMP.eXe ; rm -f VDW-TMP.eXe

# # # - 5 - M O D I F Y   F I R S T   L I N E   O F    P O S C A R 
sed -i  "1s/.*/$list1/g" $1 
sed -i  "1s/,/  /g"      $1 

