!==============================================================================
! Project: pFoF
! File: pfof_snap/src/pfof_snap.f90
! Copyright Fabrice Roy and Vincent Bouillot (2011)
! Fabrice.Roy@obspm.fr
! 
! This software is a computer program whose purpose is to detect dark matter
! haloes in cosmological simulations with a parallel Friend of Friend algorithm.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and, more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!==============================================================================

!> 
!! \mainpage pFOF Documentation
!! \section intro_sec Introduction
!!
!! pFOF is a distributed implementation of the Friends-of-Friends halo finder algorithm. <br>
!! It can be used to analyze output from the RAMSES code. <br>
!! This parallel Friend of Friend implementation is based on a sequential
!! implementation written by Edouard Audit (CEA). <br>
!! It has been written by Fabrice Roy and Vincent Bouillot -- CNRS / LUTH (Observatoire de Paris). <br>
!! mail: fabrice.roy@obspm.fr <br>
!!
!! The code has been designed, tested and optimized in collaboration with Yann Rasera. <br>
!! mail: yann.rasera@obspm.fr <br>
!! The code is described in the journal note "pFoF: a highly scalable halo-finder for large cosmological data sets", F. Roy, V. Bouillot, Y. Rasera, A&A 564, A13 (2014).
!!
!! \section install_sec Installation
!! pFOF is written in Fortran and uses MPI. <br>
!! Parallel HDF5 is required for HDF5 I/O. <br>
!! The HDF5 I/O subroutines use some Fortran 2003 features. <br>
!! A standard Makefile is provided. 
!! \subsection running Running the program
!! pFOF is an MPI application.
!! You should run it using your system MPI launcher, e.g. <BR>
!! mpirun -np 8 ./pfof 
!!
!! pFOF input parameters are read from pfof.nml text file.
!!
!! \section copyright Copyright and License
!! This license applies to etc etc...
!!
!! <BR><BR>
!!
!!

!> @file 
!! This file contains the main program of the snapshot implementation of pFoF.

Program friend

  Use modconstant
  Use mpi
  Use modmpicommons
  Use modvarcommons
  Use modvariables
  Use modhdf5
  Use modmpicom
  Use modfofpara
  Use modreadcube
  Use modwritecube
  Use modio
  Use modtiming

  Implicit none

  ! Local variables
  Integer(kind=4) :: errcode
  Integer(kind=4) :: mpierr
  
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
     Call selectreadcubetype()

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
     
     Call selectwritecubetype()

     Call writecube()
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

  If(procID==0) Then
     Call theend()
  End If

  ! Close hdf5 interface.
  Call hdf5_finalize()

  ! Finalize MPI
  Call Mpi_Finalize(mpierr)

End Program friend
