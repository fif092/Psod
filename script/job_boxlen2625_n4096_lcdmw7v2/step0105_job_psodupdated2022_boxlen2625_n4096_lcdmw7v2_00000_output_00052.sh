#!/bin/bash
#MSUB -r psod_boxlen2625_n4096_lcdmw7v2_00000
#MSUB -n 4096
#MSUB -c 2
#MSUB -M 3650
#MSUB -T 86400
#MSUB -e error_out_files/psod_boxlen2625_n4096_lcdmw7v2_00000_%I.e
#MSUB -o error_out_files/psod_boxlen2625_n4096_lcdmw7v2_00000_%I.o
#MSUB -Q normal 
#MSUB -A gen2287
#MSUB -q skylake
#MSUB -m scratch

PATHIN=/ccc/scratch/cont003/gen2287/bretonm/boxlen2625_n4096_lcdmw7v2_00000/post/ #WARNING FROM MAB FOLDER

sim=boxlen2625_n4096_lcdmw7v2_00000
POSTDIR=/ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims/$sim/post

#BE CAREFUL 200c, 500c, 2500c,DELTAVIR, ARE DIFFERENT FOR EACH REDSHIFT=>NEED TO  BE ADJUSTED BY HAND
#PLEASE ASK YANN FOR THE CODE TO COMPUTE THAT. SO YOU NEED TO LAUNCH ONE SNAP PER ONE SNAP
nummin=52
nummax=52


for ((numout=$nummin;numout<=$nummax;numout+=1)); do
    printf -v numout5 '%05d' "$numout"
    #    #WARNING ONLY 200c
    
#    #100m to 12800m WARNING REMOVED 25m and 50m AND 100M
#    for DVAL in 00200 00400 00800 01600 03200 06400 12800; do
#	echo 'Start Delta=' $DVAL 'in snap' $numout5
#	date
#	OUTDIR=$POSTDIR/sod_d${DVAL}m/output_$numout5
#	mkdir -p $OUTDIR
#	cd $OUTDIR
#	
#	ccc_mprun /ccc/cont003/home/gen2287/gen2287/proj/raygalgroupsims/code/sod/sod_profile_rdfof_parallel_hdf5 -inp $PATHIN/fof_b04000m/output_$numout5 -fil pfof_cube_snap_part_data_$sim  -nx 4096 -min 100 -sel 8 -del ${DVAL} -dl2 1#2800. -nch $numout5 -nxb 32 -nse 1000000 -sim $sim  > sod_$numout5.log 
#	
#CP info file
#	cp /ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims/boxlen2625_n4096_lcdmw7v2_00000/output_$numout5/group_00001/info_$numout5.txt .
#	
#	cd ../..
#	echo 'DONE'
#	date
#    done

#    #VIR
#    DVAL=00370.72
#    echo 'Start Delta=' $DVAL 'in snap' $numout5
#    date
#    OUTDIR=$POSTDIR/sod_deltavir/output_$numout5
#    mkdir -p $OUTDIR
#    cd $OUTDIR
#    
#    ccc_mprun /ccc/cont003/home/gen2287/gen2287/proj/raygalgroupsims/code/sod/sod_profile_rdfof_parallel_hdf5 -inp $PATHIN/fof_b04000m/output_$numout5 -fil pfof_cube_snap_part_data_$sim  -nx 4096 -min 100 -sel 8 -del ${DVAL} -dl2 1280#0. -nch $numout5 -nxb 32 -nse 1000000 -sim $sim  > sod_$numout5.log 
#    
##CP info file
#    cp /ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims/boxlen2625_n4096_lcdmw7v2_00000/output_$numout5/group_00001/info_$numout5.txt .
#    
#    cd ../..
#    echo 'DONE'
#    date

    

    #200c
    DVAL=00779.25
    echo 'Start Delta=' $DVAL 'in snap' $numout5
    date
    OUTDIR=$POSTDIR/sod_d00200c/output_$numout5
    mkdir -p $OUTDIR
    cd $OUTDIR
    
    ccc_mprun /ccc/cont003/home/gen2287/gen2287/proj/raygalgroupsims/code/sod/sod_profile_rdfof_parallel_hdf5 -inp $PATHIN/fof_b04000m/output_$numout5 -fil pfof_cube_snap_part_data_$sim  -nx 4096 -min 100 -sel 8 -del ${DVAL} -dl2 12800. -nch $numout5 -nxb 32 -nse 1000000 -sim $sim  > sod_$numout5.log 
    
#CP info file
    cp /ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims/boxlen2625_n4096_lcdmw7v2_00000/output_$numout5/group_00001/info_$numout5.txt .
    
    cd ../..
    echo 'DONE'
    date


    

#    #500c
#    DVAL=01948.12
#    echo 'Start Delta=' $DVAL 'in snap' $numout5
#    date
#    OUTDIR=$POSTDIR/sod_d00500c/output_$numout5
#    mkdir -p $OUTDIR
#    cd $OUTDIR
#    
#    ccc_mprun /ccc/cont003/home/gen2287/gen2287/proj/raygalgroupsims/code/sod/sod_profile_rdfof_parallel_hdf5 -inp $PATHIN/fof_b04000m/output_$numout5 -fil pfof_cube_snap_part_data_$sim  -nx 4096 -min 100 -sel 8 -del ${DVAL} -dl2 1280#0. -nch $numout5 -nxb 32 -nse 1000000 -sim $sim  > sod_$numout5.log 
#    
##CP info file
#    cp /ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims/boxlen2625_n4096_lcdmw7v2_00000/output_$numout5/group_00001/info_$numout5.txt .
#    
#    cd ../..
#    echo 'DONE'
#    date




#    #2500c
#    DVAL=09740.60
#    echo 'Start Delta=' $DVAL 'in snap' $numout5
#    date
#    OUTDIR=$POSTDIR/sod_d02500c/output_$numout5
#    mkdir -p $OUTDIR
#    cd $OUTDIR
#    
#    ccc_mprun /ccc/cont003/home/gen2287/gen2287/proj/raygalgroupsims/code/sod/sod_profile_rdfof_parallel_hdf5 -inp $PATHIN/fof_b04000m/output_$numout5 -fil pfof_cube_snap_part_data_$sim  -nx 4096 -min 100 -sel 8 -del ${DVAL} -dl2 1280#0. -nch $numout5 -nxb 32 -nse 1000000 -sim $sim  > sod_$numout5.log 
#    
##CP info file
#    cp /ccc/scratch/cont003/gen2287/gen2287/raygalgroupsims/boxlen2625_n4096_lcdmw7v2_00000/output_$numout5/group_00001/info_$numout5.txt .
#    
#    cd ../..
#    echo 'DONE'
#    date
    

    
done



echo 'END'

