#!/bin/bash
#MSUB -r haloanalyzersod_boxlen2625_n4096_lcdmw7v2_00000
#MSUB -n 4096
#MSUB -c 1
#MSUB -M 3650
#MSUB -T 86400
#MSUB -e error_out_files/haloanalyzersod_boxlen2625_n4096_lcdmw7v2_00000_%I.e
#MSUB -o error_out_files/haloanalyzersod_boxlen2625_n4096_lcdmw7v2_00000_%I.o
#MSUB -Q normal 
#MSUB -A gen2287
#MSUB -q skylake
#MSUB -m scratch

set -x

echo 'BEGIN'

sim=boxlen2625_n4096_lcdmw7v2_00000
PROJDIR=/ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims
POSTDIR=${PROJDIR}/$sim/post
DATADIR=${PROJDIR}/data/

nummin=52
nummax=52

nfiles=4096 #SHOULD MATCH MSUB -n NUMBER

for ((numout=$nummin;numout<=$nummax;numout+=1)); do
    printf -v numout5 '%05d' "$numout"

     #100m to 12800m (ONLY) WARNING REMOVED 25m and 50m and 100m (after 41)
    #for DVALFULL in  00200m 00400m 00800m 01600m 03200m 06400m 12800m eltavir 00200c 00500c 02500c; do #WARNING ONLY 200c
    for DVALFULL in 00200c; do
	echo 'Start Delta=' $DVALFULL 'in snap' $numout5
	date
	OUTDIR=$POSTDIR/sod_d${DVALFULL}/output_$numout5
	mkdir -p $OUTDIR
	cd $OUTDIR

	#CREATEDIR	
	rm filelist.nml
   
	cat > filelist.nml << EOF
&Size
filelistsize=${nfiles} /

&List
EOF
	
	for ((i=0;i<${nfiles}-1;i+=1)); do
	    N=$(printf %5.5i $i) #output to analyze
	    j=`expr $i + 1`      #fortran indices start from 1
	    echo "filelist($j) = '${OUTDIR}/psod_halo_snap_part_data_${sim}_${N}.h5'" >> filelist.nml
	done
	i=$(echo $nfiles '-' 1 |bc -l) 
	N=$(printf %5.5i $i) #output to analyze                                                                                                              
	j=`expr $i + 1`      #fortran indices start from 1                                                                                                         
	echo "filelist($j) = '${OUTDIR}/psod_halo_snap_part_data_${sim}_${N}.h5' / " >> filelist.nml
	
	
	
	

	#JOB ITSELF
	
        #parameters
	if [ $DVALFULL = '00050m' ]; then
	    DVAL='00050'
	fi
	if [ $DVALFULL = '00100m' ]; then
	    DVAL='00100'
	fi
	if [ $DVALFULL = '00200m' ]; then
	    DVAL='00200'
	fi
	if [ $DVALFULL = '00400m' ]; then
	    DVAL='00400'
	fi
	if [ $DVALFULL = '00800m' ]; then
	    DVAL='00800'
	fi
	if [ $DVALFULL = '01600m' ]; then
	    DVAL='01600'
	fi
	if [ $DVALFULL = '03200m' ]; then
	    DVAL='03200'
	fi
	if [ $DVALFULL = '06400m' ]; then
	    DVAL='06400'
	fi
	if [ $DVALFULL = '12800m' ]; then
	    DVAL='12800'
	fi


	#SPECIFIC TO A GIVEN REDSHIFT
	if [ $DVALFULL = 'eltavir' ]; then
	    DVAL='00370.72'
	fi
	if [ $DVALFULL = '00200c' ]; then
	    DVAL='00779.25'
	fi
	if [ $DVALFULL = '00500c' ]; then
	    DVAL='01948.12'
	fi
	if [ $DVALFULL = '02500c' ]; then
	    DVAL='09740.60'
	fi
	

	delta=$DVAL    #SOD OVERDENSITY
	data=$DATADIR
	
	
        #get simu name, output number and box length                      
	echo $sim
	cosmo=$(echo ${sim} | awk -F'_' '{ print $3 }')
	boxlen=$(echo ${sim} | awk -F'_' '{ print $1 }')
	boxlength=$(echo ${boxlen} | awk -F'len' '{ print $2 }')
	
#create parameter file
	cat > haloanalyzersod_$sim.param <<EOF
${boxlength}
$delta
$data/ramses_input_${cosmo}.dat
EOF
	
#RUN
	ccc_mprun /ccc/cont003/home/gen2287/gen2287/proj/raygalgroupsims/code/haloanalyzer_new/src/haloanalyzersod haloanalyzersod_${sim}.param >haloanalyzersod_${sim}.log



	cd ../..
	echo 'DONE'
	date

    done
    
done

echo 'END'


