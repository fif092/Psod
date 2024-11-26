#!/bin/bash

usage() {
    echo "Create hdf5 cone files from Ramses output_ncoarse files."
    echo "mandatory arguments:"
    echo "-d /PATH/TO/RAMSES/OUTPUT_NCOARSE_DIRS"
    echo "-f NameOfTheConeFiles"
    echo "-o OutputFileName"
    echo "-c CubeSize"
    echo "optional arguments:"
    echo "-e ExecName (default: ./conecreator)"
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
EXEC=./conecreator

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
	    thetay=`grep thetay $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	    THETAYCONE=$(echo $thetay | awk -F= '{print $2}')
	    thetaz=`grep thetaz $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	    THETAZCONE=$(echo $thetaz | awk -F= '{print $2}')
	    nstride=`grep nstride $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	    NSTRIDERAMSES=$(echo $nstride | awk -F= '{print $2}')
	    fs=`grep isfullsky $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	    FSSTRING=$(echo $fs | awk -F= '{print $2}')
	    if [ $FSSTRING == '0' ]
	    then
		FSSTRING='.false.'
	    else
		FSSTRING='.true.'
	    fi
	    echo "first info file found: "$fileinfo
	    echo "height of the cone or radius of the sphere (dmax): "$RMAXCONE
	    echo "fullsky :"$FSSTRING
	    echo " "
	fi
	npart=`grep nglobalcell $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	NPARTSHELL=$(echo $npart | awk -F= '{print $2}')
	ncpu=`grep ncpu $DIR$dirncoarse'/'$fileinfo | sed 's/ //g'`
	NFILEMAX=$(echo $ncpu | awk -F= '{print $2}')
	OUTPUTNAME=$OUT'_'$ncoarse'.h5'

	if [ -f infocone.nml ] 
	then
	    rm -f infocone.nml
	fi 
	
	cat >> infocone.nml <<EOF
&infocone
rmax   = $RMAXCONE
thetay = $THETAYCONE
thetaz = $THETAZCONE
npart  = $NPARTSHELL
fullsky = $FSSTRING
/

&infofile
dirname  = '$DIR$dirncoarse'
filename = '${FILE}_${ncoarse}_proc'
filenb   = $nfile
ffile    = 1
lfile    = $NFILEMAX
nstride  = $NSTRIDERAMSES
/

&infooutput
shellname = '$OUTPUTNAME'
cubesize = $CUBE
/
EOF
	echo "conecreator run in directory "$dirncoarse
	mpirun -np ${NCPU} ${EXEC}
	echo " "
    fi
done

LID=$ncoarse

echo " "
echo "conemapper run using shell files "$OUT
echo "first shell file: "$FID
echo "last shell file: "$LID
echo " "
if [ -f pfof_cone.nml ] 
then
    rm -f pfof_cone.nml
fi 

cat >> pfof_cone.nml <<EOF
&shellparameters
shell_dir = '$PWD/'
shell_filename ='$OUT'
shell_fid = $FID
shell_lid = $LID
/

&fofparameters
perco = 0.20
Mmin = 100
Mmax = 100000000
nres = 128
/
 
&outputparameters
output_root = 'test'
/

EOF
./conemapper
