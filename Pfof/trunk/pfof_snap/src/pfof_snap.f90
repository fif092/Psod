! Copyright (c) 2007-2016 CNRS, Fabrice Roy, Vincent Bouillot
! Author: Fabrice Roy (LUTH/CNRS/PSL), fabrice.roy@obspm.fr
! Vincent Bouillot (LUTH/CNRS/PSL)
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

!> 
!! \mainpage pFOF Documentation
!! \section intro_sec Introduction
!!
!! pFOF is a software package based on a distributed implementation of the Friends-of-Friends halo finder algorithm. <br>
!! It can be used to analyze output from the RAMSES code. <br>
!! This parallel Friend of Friend implementation is based on a sequential
!! implementation written by Edouard Audit (CEA). <br>
!! It has been written by Fabrice Roy and Vincent Bouillot -- CNRS / LUTH (Observatoire de Paris). <br>
!! mail: fabrice.roy@obspm.fr <br>
!! The tools have been written by Fabrice Roy. <br>
!!
!! The package has been designed, tested and optimized in collaboration with Yann Rasera. <br>
!! mail: yann.rasera@obspm.fr <br>
!! The first version of the code is described in the journal note "pFoF: a highly scalable halo-finder for large cosmological data sets", F. Roy, V. Bouillot, Y. Rasera, A&A 564, A13 (2014). <br>
!!
!! \section content_sec Content of the package
!! pFOF contains several softwares:
!! - pfof_snap: 'snapshot' version, used to detect dark matter halos in Ramses snapshots;
!! - pfof_cone: 'cone' version, used to detect dark matter halos in Ramses lightcones;
!! - conepartcreator: creates hdf5 containing shells of particles from Ramses output_ncoarse files;
!! - conegravcreator: creates hdf5 containing shells of cells from Ramses output_ncoarse files;
!! - conemapper: creates a mapping files from hdf5 shells to prepare a pfof_cone analyzis;
!! - regionextractor: extract the particles from cubes files (created by pfof_snap) given a position and the size of the region to extract;
!! - haloanalyzer; skeleton of a tool that can be used to further analysis of the halos detected by pfof_snap or pfof_cone.
!! 
!! \section install_sec Installation
!! pFOF is written in Fortran 2003 and uses MPI. <br>
!! Parallel HDF5 is required for HDF5 I/O. <br>
!! The HDF5 I/O subroutines use some Fortran 2003 features. <br>
!! HDF5 should be configured with the flags --enable-parallel and --enable-fortran (and --enable-fortran2003 for versions < 1.10.0).<br>
!! A standard Makefile is provided for each software in the package. 
!! 
!! \subsection running Running the program
!! Running pfof_snap to analyze Ramses snapshots is straightforward:<br>
!! mpirun -np 8 ./pfof<br>
!! Note that you should use a cubic power as the number of processes (i.e. 8, 27, 64, 125, 216, etc.).<br>
!!
!! pfof_snap input parameters are read from pfof_snap.nml text file.<br>
!!
!! To run a lightcone analysis you should first create the shell files with conepartcreator, the create the mapping with conemapper, then copy the mapping that requires an adequate number of processes to
!! procmap.h5, then run pfof_cone with the required number of processes.
!!
!! 
!!
!! \subsection test_sec Test
!!
!! You can download some data set to test your installation.<br>
!! You can run every software in the package with this dataset.<br>
!! You just have to modify the pathes in the <br>
!!
!! \section copyright Copyright and License
!! This license applies to etc etc...
!!
!! <BR><BR>
!!
!!

!> @file pfof_snap.f90
!!
!! @author Fabrice Roy
!! @author Vincent Bouillot
!!
!! This file contains the main program.
!!
Program friend

  Use modhdf5
  Use modvariables
  Use modparameters
  Use modmpicom
  Use modfofpara
  Use modreadcube
  Use modwritecube
  Use modio
  Use modtiming

  Implicit none

  ! Local variables
  Integer :: errcode

  !----------------------------------------------------

  ! Initialization of MPI
  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)

  ! Print "title"
  If(procID==0) Then
     Call title()
  End If
  
  ! Read parameters
  Call readparameters()
  
  ! Print parameters on screen
  If(procID==0) Then
     Call print_screen_parameters()
  End If

  ! Creates the different communicators used for the analysis and the // I/O
  Call setcommunicators()
  
  ! Initialize HDF5 interface
  Call hdf5_init()


  ! Particles read from previously created cube files
  If( param%do_read_from_cube ) Then
     Call selectcubetype()

     If(param%do_timings) Then
        timeInt = Mpi_Wtime()
     End If
     Call readcube()
     If(param%do_timings) Then
        tRead = Mpi_Wtime() - timeInt
        tReadfile = tRead
     End If
  Else If (param%code_index.eq.'RA2' .or. param%code_index.eq.'RA3') Then
     ! code_index determines the type of cosmological simulation to be analyzed:
     ! RA2 = RAMSES v2
     ! RA3 = RAMSES v3
     ! Read the output of the cosmological simulation
     Call RAMSES_lecture()
  Else
     Print *,'Proc',procID, '  TYPE',param%code_index
     If(procID==0) Print *,'** Wrong file type. Possibilities are RA2 or RA3. **'
     errcode = 1
     Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
     Stop
  End If
  
  ! Print on screen and in log file 
  If(procID==0) Then
     Write(* ,*)   'nz = ',nres,'ngrid = ',ngrid
     Write(Ulog,*) 'nz = ',nres,'ngrid = ',ngrid
  End If


  ! Write cube of particles files if requested
  If(param%do_write_cube) Then
     If(procID==0) Print *,      'Write particles distributed in a cartesian grid'
     If(procID==0) Write(Ulog,*) 'Write particles distributed in a cartesian grid'
     
     If(param%gatherwrite_factor==1) Then
        If(param%do_sort_cube) Then
           Call h5writesortedcube()
        Else
           Call h5writecube()
        End If
     Else
        If(param%do_sort_cube) Then
           Call mpih5writesortedcube()
        Else
           Call mpih5writecube()
        End If
     End If
  End If


  If(procID==0) Then
     Print *,' '
     Print *,'Friends of Friends halo detection'
     Print *,' '
     Write(Ulog,*) 'Friends of Friends halo detection'
  End If


  ! Parallele Friends of Friends if requested
  If(param%do_fof) Then
     Call fofpara()
  Else
     tFoF     = 0.0
     tFoFinit = 0.0
     tFoFloc  = 0.0
     tRaccord = 0.0
     tObs     = 0.0
     tOut     = 0.0
  End If

  
  ! Output of timings if requested
  If(param%do_timings .and. procID==0) Then
     Call writetimings()
  End If
  

  ! Close hdf5 interface.
  Call hdf5_finalize()

  If(procID==0) Then
     Call theend()
  End If

  ! Finalize MPI
  Call Mpi_Finalize(mpierr)

End Program friend
