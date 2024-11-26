!> 
!! \mainpage pFOF Documentation
!! \section intro_sec Introduction
!!
!! pFOF is a distributed implementation of the Friends-of-Friends halo finder algorithm. <br>
!! It can be used to analyze output from the RAMSES code. <br>
!! This parallel Friend of Friend implementation is based on a sequential
!! implementation written by Edouard Audit (CEA). <br>
!! It has been written by Fabrice Roy -- CNRS / LUTH (Observatoire de Paris). <br>
!! mail: fabrice.roy@obspm.fr <br>
!! An alternative version that can also analyze lightcones simulated with the RAMSES code has been developped by Vincent Bouillot. <br>
!! mail: vincent.bouillot@obspm.fr <br>
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
!! This file contains the main program.
!!
Program friend

#ifdef WITHHDF5
  Use modhdf5
#endif
  Use modconst
  Use modvariable
  Use modparam
  Use modmpicom
  Use modfofpara
  Use modio
  Use modtiming

  Implicit none

  ! Local variables
  Integer :: errcode
  Character(len=15) :: codeversion
  Integer :: ierr

  !----------------------------------------------------

  ! Initialization of MPI
  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)

  If(procID==0) Then
     Call title()
  End If
  
  Call readparameters()
  
  If(procID==0) Then
     Call writeparameters()
  End If

  Call setcommunicators()
  
#ifdef WITHHDF5
  If(usehdf5) Then
     ! Open hdf5 interface
     Call hdf5_init()

     Open(Unit=10, File='pfof.version', status='old', Iostat=ierr)
     If(ierr/=0) Then
        Print *, 'version of pfof not defined: file not found'
        codeversion = 'undefined'
     Else
        Read(10,*,Iostat=ierr) codeversion
        If(ierr/=0) Then
           Print *, 'version of pfof not defined: file empty'
           codeversion='undefined'
        End If
        Close(10)
     End If
     origin = 'Created with pfof version '//codeversion
  End If
#endif
  
  ! code_index determines the type of cosmological simulation to be analyzed:
  ! RA2 = RAMSES v2
  ! RA3 = RAMSES v3
  ! Read the output of the cosmological simulation

  If( readfromcube ) Then
#ifdef WITHHDF5
     If(usehdf5) Then
        If(gatherread==1) Then
           Call h5readcube()
        Else
           Call mpih5readcube()
        End If
     Else
        Print *,'Read from cube is only available with hdf5 i/o file format. Parameter usehdf5 should be .true. '
        errcode = 2
        Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
     End If
#else
     Print *,'Read from cube is only available with hdf5 i/o file format.'
     errcode = 2
     Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
#endif
  Else If (code_index.eq.'RA2' .or. code_index.eq.'RA3') Then
     Call RAMSES_lecture()
  Else
     If(procID==0) Print *,'** Wrong file type. Possibilities are RA2 or RA3. **'
     errcode = 2
     Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
     Stop
  End If
  
  ! Print on screen and in log file 
  If(procID==0) Then
     Write(* ,*)   'nz = ',nres,'ngrid = ',ngrid
     Write(Ulog,*) 'nz = ',nres,'ngrid = ',ngrid
  End If




  ! Write cube of particles files if requested
  If(outcube) Then
     If(procID==0) Print *,      'Write particles distributed in a cartesian grid'
     If(procID==0) Write(Ulog,*) 'Write particles distributed in a cartesian grid'
     
#ifdef WITHHDF5
     If(usehdf5) Then
        If(gatherwrite==1) Then
           Call h5writecube()
        Else
           Call mpih5writecube()
        End If
     Else
        Call outputcube()
     End If
#else
     Call outputcube()
#endif

  End If


  If(procID==0) Then
     Print *,' '
     Print *,'Friends of Friends halo detection'
     Print *,' '
     Write(Ulog,*) 'Friends of Friends halo detection'
  End If


  ! Parallele Friends of Friends if requested
  If(dofof) Then
     Call parafof()
  Else
     tFoF     = 0.0
     tFoFinit = 0.0
     tFoFloc  = 0.0
     tRaccord = 0.0
     tObs     = 0.0
     tOut     = 0.0
  End If

  
  ! Output of timings if requested
  If(dotimings .and. procID==0) Then
     Call writetimings()
  End If
  

#ifdef WITHHDF5
  If(usehdf5) Then
     ! Close hdf5 interface.
     Call hdf5_finalize()
  End If
#endif

  ! Finalize MPI
  Call Mpi_Finalize(mpierr)

End Program friend
