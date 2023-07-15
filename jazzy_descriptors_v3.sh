#!/bin/bash

#SCRIPT TO OBTAIN HYDRATION ENERGIES WITH JAZZY

#table must contain:
#Column 1 -> id
#Column 2 -> SMILES
#Column 3 -> logP data

########
#INPUTS#
########

id=$1
table=$2

var0=`cut -d";" -f 1 $table`
ids=($var0)

var=`cut -d";" -f 2 $table`
smiles=($var)

l0=`echo "${#smiles[@]}"`
l=`echo $[$l0 - 1]`

for i in $(eval echo "{1..$l}")
do
	
id1=`echo ${ids[$i]}`

CHds=`jazzy vec --opt MMFF94 ${smiles[$i]} | awk '{print $2}' | cut -d"," -f 1`
XHds=`jazzy vec --opt MMFF94 ${smiles[$i]} | awk '{print $4}' | cut -d"," -f 1`
HBAs=`jazzy vec --opt MMFF94 ${smiles[$i]} | awk '{print $6}' | cut -d"," -f 1`

HydA0=`jazzy vec --opt MMFF94 ${smiles[$i]} | awk '{print $8}' | cut -d"," -f 1`
HydA=`echo "$HydA0/4.184" | bc -l`
HydP0=`jazzy vec --opt MMFF94 ${smiles[$i]} | awk '{print $10}' | cut -d"," -f 1`
HydP=`echo "$HydP0/4.184" | bc -l`
Hyd0=`jazzy vec --opt MMFF94 ${smiles[$i]} | awk '{print $12}' | cut -d"," -f 1`
Hyd=`echo "$Hyd0/4.184" | bc -l`

echo "ID,smiles,CH_strength,XH_strength,HBA_strength,Hyd_Apolar,Hyd_Polar,Hyd">Title.csv
echo "$id1,$id2,$CHds,$XHds,$HBAs,$HydA,$HydP,$Hyd" >> Data.csv
cat Title.csv Data.csv > jazzy_descriptors.$id.csv

done 

rm Title.csv
rm Data.csv
