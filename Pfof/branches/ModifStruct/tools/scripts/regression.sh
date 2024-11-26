#!/bin/bash

REGRESSIONDIR=$HOME/Cosmo/test/TestRegression
TRUNKDIR=$HOME/Cosmo/Pfof/trunk
RELEASEDIR=$HOME/Cosmo/Pfof/tags/Release-1.0

# $1: name of the code to test: if 'all' then each test is called
main() {
    if (($#==0))
    then local testname='all'
    else local testname=$1
    fi
    
    echo " "
    echo "---------------------------------------"
    echo " Regression tests for the pfof package "
    echo "---------------------------------------"
    echo " "

    case "$testname" in

	all)  echo "Performing all regression tests"
	    testall
	    ;;
	pfofhdf5)  echo "Performing pfofhdf5 regression test"
	    testpfofhdf5
	    ;;
	pfofcone)  echo "Performing pfofcone regression test"
	    testpfofcone
	    ;;
	conecreator) echo "Performing conecreator regression test"
	    testconecreator
	    ;;
	conemapper) echo "Performing conemapper regression test"
	    testconemapper
	    ;;
	haloanalyzer) echo "Performing haloanalyzer regression test"
	    testhaloanalyzer
	    ;;
	*) usage
	    ;;
	
esac
    
}

usage() {
    echo "Regression test for the pfof package"
    echo "Choose a test from all, conecreator, conemapper, halohanalyzer, pfofcone, pfofhdf5"
    echo "Usage: $0 testname"
    echo " "
}

testall() {
    testpfofhdf5
    testpfofcone
    testconecreator
    testconemapper
    testhaloanalyzer
}

diffpfof() {

    DATASETHALOOK=true
    DATASETCUBEOK=true
    DATASETMASSOK=true
    DATASETMPICOK=true

    for file in `ls regtrunk*.h5`
    do
	FILETRUNK=$file
	FILERELEASE=$(echo $file | sed -e 's/trunk/release/g')
	#echo "Trunk: $FILETRUNK, release: $FILERELEASE"
	H5TYPE=$(echo $file | awk -F '_' '{ print $2 }' )
	case "$H5TYPE" in
	    halo) DIFFDATA=`h5diff $FILETRUNK $FILERELEASE | grep 'dataset'`
		if [ $? -eq 1 ] && [ "$DATASETHALOOK"=true ]
		then
		    DATASETHALOOK=true
		else
		    DATASETHALOOK=false
		fi
		;;
	    cube) DIFFDATA=`h5diff $FILETRUNK $FILERELEASE | grep 'dataset'`
		if [ $? -eq 1 ] && [ "$DATASETCUBEOK"=true ]
		then
		    DATASETCUBEOK=true
		else
		    DATASETCUBEOK=false
		fi
		;;
	    halomass.h5) DIFFDATA=`h5diff $FILETRUNK $FILERELEASE | grep 'dataset'`
		if [ $? -eq 1 ] && [ "$DATASETMASSOK"=true ]
		then
		    DATASETMASSOK=true
		else
		    DATASETMASSOK=false
		fi
		;;
	    mpicube) DIFFDATA=`h5diff $FILETRUNK $FILERELEASE | grep 'dataset'`
		if [ $? -eq 1 ] && [ "$DATASETMPICOK"=true ]
		then
		    DATASETMPICOK=true
		else
		    DATASETMPICOK=false
		fi
		;;

	esac
    done

    if [ "$DATASETHALOOK"=true ] 
    then
	echo "Halo files: OK"
    else
	echo "Halo files: some datasets differ"
    fi
    if [ "$DATASETCUBEOK"=true ] 
    then
	echo "Cube files: OK"
    else
	echo "Cube files: some datasets differ"
    fi
    if [ "$DATASETMASSOK"=true ] 
    then
	echo "Mass files: OK"
    else
	echo "Mass files: some datasets differ"
    fi
    
    echo " "

}

testpfofhdf5() {
    echo "Warning: pfofhdf5 regression test not fully implemented"
    echo " "

    cd $TRUNKDIR
    cd pfof_hdf5/src
    TRUNKVERSION=`svnversion`
    echo "SVN version of trunk: $TRUNKVERSION"
    echo "Build pfofhdf5"
    echo " "
    make
    echo " "
    
    cp pfofhdf5 $REGRESSIONDIR/pfofhdf5trunk
    cp ../pfof.nml $REGRESSIONDIR/pfoftrunk.nml
    
    cd $RELEASEDIR
    cd pfof_hdf5/src
    RELEASEVERSION=`svnversion`
    echo "SVN version of release: $RELEASEVERSION"
    echo "Build pfofhdf5"
    echo " "
    make 
    echo " "

    cp pfofhdf5 $REGRESSIONDIR/pfofhdf5release
    cp ../pfof.nml $REGRESSIONDIR/pfofrelease.nml

    cd $REGRESSIONDIR

    ###
    echo "First test: read from Ramses, HDF5 outputs, b=0.2, Mmin=100" 
    OLDROOT=$(cat pfoftrunk.nml | grep root)
    OUTPUTROOT=$(echo $OLDROOT | awk -F '=' '{ print $1 }' | sed -e 's/\x27//g')
    OUTPUTROOT="$OUTPUTROOT = 'regtrunk'"
    sed -i.bak s/"$OLDROOT"/"$OUTPUTROOT"/g pfoftrunk.nml 
    cp pfoftrunk.nml pfof.nml

    echo "time mpirun -np 8 ./pfofhdf5trunk"
    time mpirun -np 8 ./pfofhdf5trunk > trunk1.log
    echo " "

    OLDROOT=$(cat pfofrelease.nml | grep root)
    OUTPUTROOT=$(echo $OLDROOT | awk -F '=' '{ print $1 }' | sed -e 's/\x27//g')
    OUTPUTROOT="$OUTPUTROOT = 'regrelease'"
    sed -i.bak s/"$OLDROOT"/"$OUTPUTROOT"/g pfofrelease.nml
    cp pfofrelease.nml pfof.nml
    
    echo "time mpirun -np 8 ./pfofhdf5release"
    time mpirun -np 8 ./pfofhdf5release > release1.log
    echo " "
    
    diffpfof

    ###
    echo "Second test: read from HDF5 cubes, HDF5 gathered outputs, b=0.2, Mmin=100" 

    READCUBE=$(cat pfoftrunk.nml | grep readfromcube | sed -e 's/\///g')
    sed -i.bak s/"$READCUBE"/"readfromcube = .true."/g pfoftrunk.nml 
    GATHERWRITE=$(cat pfoftrunk.nml | grep gatherwrite )
    sed -i.bak s/"$GATHERWRITE"/"gatherwrite = 2"/g pfoftrunk.nml
    cp pfoftrunk.nml pfof.nml

    echo "time mpirun -np 8 ./pfofhdf5trunk"
    time mpirun -np 8 ./pfofhdf5trunk > trunk2.log
    echo " "
 
    READCUBE=$(cat pfofrelease.nml | grep readfromcube | sed -e 's/\///g')
    sed -i.bak s/"$READCUBE"/"readfromcube = .true."/g pfofrelease.nml
    GATHERWRITE=$(cat pfofrelease.nml | grep gatherwrite )
    sed -i.bak s/"$GATHERWRITE"/"gatherwrite = 2"/g pfofrelease.nml
    cp pfofrelease.nml pfof.nml

    echo "time mpirun -np 8 ./pfofhdf5release"
    time mpirun -np 8 ./pfofhdf5release > release2.log
    echo " "

    diffpfof

    ###
    echo "Third test: read from HDF5 gathered cubes, HDF5 outputs, b=0.3, Mmin=200"
}

testpfofcone() {
    echo "pfofcone regression test not implemented"
    echo " "
}

testconecreator() {
    echo "conecreator regression test not implemented"
    echo " "
}

testconemapper() {
    echo "conemapper regression test not implemented"
    echo " "
}

testhaloanalyzer() {
    echo "haloanalyzer regression test not implemented"
    echo " "
}

main "$@"

#make $(EXE)
#	mkdir -p $(REGDIR)
#	cp $(EXE) $(REGDIR)
#	cp ../pfof.nml $(REGDIR)
#	cd $(REGDIR) && pwd && NAME=$(cat pfof.nml | grep root | awk -F '=' '{ print $2 }' | sed -e 's/\x27//g')
#	echo $(NAME)
#	#mpirun -np 8 ./pfofhdf5
#	cd $(RELEASEDIR)/pfof_hdf5
#	make $(EXE)
#	cp $(EXE) $(REGDIR)/pfofhdf5_release
#	cp ../pfof.nml $(REGDIR)
#	cd $(REGDIR)
	#mpirun -np 8 ./pfofhdf5_release


