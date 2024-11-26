#!/bin/bash

usage() {
    echo "Create hdf5 particles cone files from Ramses output_ncoarse files."
    echo "mandatory arguments:"
    echo "-d /PATH/TO/RAMSES/OUTPUT_NCOARSE_DIRS"
    echo "-f NameOfTheConeFiles"
    echo "-o OutputFileName"
    echo "-c CubeSize"
    echo "optional arguments:"
    echo "-e ExecName (default: ./conepartcreator)"
    echo "-n Ncpu (default: 2)"
}

echo " "
echo "----------------------------------------------------------------------"
echo "Creation of the hdf5 shell files and the process mapping for pfof_cone"
echo "----------------------------------------------------------------------"
echo " "

if (($#<8))
    then usage
    exit 1
fi

# Default value for optional arguments
NCPU=2
EXEC=./conepartcreator

while getopts ":h:d:f:o:c:e:n:" Option
do
    case $Option in
        h ) usage;;
        d ) DIR=$OPTARG;;
        f ) FILE=$OPTARG;;
        o ) OUT=$OPTARG;;
        c ) CUBE=$OPTARG;;
	e ) EXEC=$OPTARG;;
	n ) NCPU=$OPTARG;;
        * ) echo "Unimplemented option chosen. Try $0 -h";
            exit 1;;   # DEFAULT
    esac
done

echo "Options used"
echo "location of the Ramses output_ncoarse directories: "$DIR
echo "name of the Ramses cone files: "$FILE
echo "name of the shell output files: "$OUT
echo "size of the discretization cubes used for the shell files: "$CUBE
echo " "

RMAXCONE='0'
FID='0'
LID='0'

for dirfullpath in `ls -d $DIR/output_ncoarse*`;do
    #echo $dirfullpath
    dirncoarse=$(echo $dirfullpath | awk -F/ '{print $NF}')
    #echo $dirncoarse
    ncoarse=$(echo $dirncoarse | awk -F_ '{print $3}')
    #echo $ncoarse
    nfile=`ls -1 $dirfullpath/$FILE*.dat 2>/dev/null | wc -l`
    #echo $nfile
    if [ $nfile -ne 0 ] 
    then
	if [ $FID == '0' ] 
	    then
	    FID=$ncoarse
	fi
	fileinfo='info_'$FILE'_'$ncoarse'.txt'
	if [ $RMAXCONE == '0' ]
	then
	    dmax=`grep dmax $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	    RMAXCONE=$(echo $dmax | awk -F= '{print $2}')
	    echo "first info file found: "$fileinfo
	    echo "height of the cone or radius of the sphere (dmax): "$RMAXCONE
	    echo " "
	fi
	ncpu=`grep ncpu $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	NFILEMAX=$(echo $ncpu | awk -F= '{print $2}')

	echo "ncpu : "$NFILEMAX
	if [ -f infoconepart.nml ] 
	then
	    rm -f infoconepart.nml
	fi 
	
	cat >> infoconepart.nml <<EOF
&input_parameters
input_path = '$DIR$dirncoarse'
cone_input_file = '${FILE}_${ncoarse}_proc'
info_ramses_input_file = 'info_ncoarse_${ncoarse}.txt'
info_cone_input_file = 'info_${FILE}_${ncoarse}.txt'
nfile = $nfile
first_file = 1
last_file = $NFILEMAX
do_read_ramses_part_id = .true. 
do_read_potential = .true. 
do_read_gravitational_field = .true. 
cone_max_radius = $RMAXCONE
/

&output_parameters
simulation_name = '$OUT'
cube_size = $CUBE
/
EOF
	cp infoconepart.nml infoconepart_${ncoarse}.nml
	echo "conecreator run in directory "$dirncoarse
	mpirun -np ${NCPU} ${EXEC}
	echo " "
    fi
done

exit 0

