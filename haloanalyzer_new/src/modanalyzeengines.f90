!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================
Module modanalyzeengines

  Use modvariables
  Use modfunctions
  Use modhdf5
  Implicit None

  Type ListOfVar
     Logical :: pos = .true.
     Logical :: vel = .true.
     Logical :: pot = .true.
     Logical :: id  = .true.
  End type ListOfVar
  
  Abstract interface 
     Subroutine generic_function()
     End Subroutine generic_function
  End interface

Contains

  !=======================================================================
  !> Analyze a list of halos defined in the file halolistfilename
  !! Apply the function generic_function to each halo individually
  Subroutine analyzelist(halolistfilename, myfunction, inputvarlist)

    Use modiocommons

    Character(len=400), intent(in) :: halolistfilename
    Procedure(generic_function) :: myfunction
    Type(ListOfVar), intent(in), optional :: inputvarlist
    Type(ListOfVar)     :: varlist
    Character(len=391)  :: basename        ! name of the file without the number of the file and the extension
    Character(len=400)  :: filename        ! name of the file to open
    Integer :: i, ih, ifile                ! loop indices
    Integer :: nfile                       ! number of halo files written by pFoF 

    Integer :: fh, lh
    Integer :: halolistsize                ! number of halo in the list of halo read from the parameter file
    Integer :: myhalolistsize

    Integer(kind=PRI), dimension(:), allocatable :: halolist    ! list of the ID of the halo that we must read
    Namelist / Size / halolistsize    ! namelist
    Namelist / File / basename        
    Namelist / List / halolist        ! namelist

    If(Present(inputvarlist)) Then
       varlist = inputvarlist
    End If

#ifdef WITHMPI
    Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
    Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
    
    If(procID == 0) Then
#endif

       ! We read the list of halo that we want to read and analyze: their number and ID
       Open(Unit=10,file=trim(halolistfilename),status='old')
       Read(10,nml=Size)
       Read(10,nml=File)
       Allocate(halolist(halolistsize+1))
       Read(10,nml=List)
       Close(10)

       ! We add a fake haloID at the end of the halo list
       halolist(halolistsize+1) = 0

       fh = 1
       lh = halolistsize
       myhalolistsize = halolistsize

#ifdef WITHMPI
    End If
    
    Call Mpi_Bcast(basename, len(basename), Mpi_Char, 0, Mpi_Comm_World, mpierr)
    Call Mpi_Bcast(halolistsize, 1, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
    If(procID /= 0) Then
       Allocate(halolist(halolistsize+1))
    End If

    myhalolistsize = halolistsize / procNB
    fh = procID*myhalolistsize + 1
    If( procID < mod(halolistsize, procNB) ) Then
       myhalolistsize = myhalolistsize + 1
       fh = fh + procID
    Else
       fh = fh + mod(halolistsize,procNB)
    End If

    lh = fh + myhalolistsize - 1

    Call Mpi_Bcast(halolist, halolistsize+1, MPI_PRI, 0, Mpi_Comm_World, mpierr)
    
#endif


!    basename = 'test_halo'
    filenbchar = '_00000.h5'
    nfile = 1

    nbhaloanalyzed = myhalolistsize

    ! We assume that there is at least 1 file.
    ! We read the total number of files written by pFoF in this first file.
    ! We read what is to be read in this file, and move to the next and repeat until the last file is read. 
    ifile = 0
    ih = fh 
    Do While(ifile < nfile  .And. ih <= lh )

       Write(filenbchar(2:6),'(I5.5)') ifile
       filename = trim(basename)//trim(filenbchar)
       Call hdf5_open_file(filename, file_id)
       name = 'nfile'
       Call hdf5_read_attr(file_id, name, nfile)
       name = 'nhalo_file'
       Call hdf5_read_attr(file_id, name, haloNB)
       name = 'identity_halo'
       Allocate(haloID(haloNB))
       Call hdf5_read_data(file_id, name, haloNB, haloID)

       ! We assume that the halo list is sorted 
       ! If the ID of the 1st halo we are looking for is greater than the ID of the last halo in the file
       ! we move to the next file
       If( halolist(ih) > haloID(haloNB) ) Then
          ifile = ifile + 1
          Deallocate(haloID)
          Cycle
       Else 
          Do While( halolist(ih) <= haloID(haloNB) .And. ih <= lh  )

             currenthaloID = halolist(ih)
             currenthalo = ih - fh + 1

             groupname = 'halo_0000000000000000000'
             Write(groupname(6:24),'(I19.19)') currenthaloID
             Call hdf5_open_group(gr_id, groupname, gr_halo_id)
             name = 'npart_halo'
             Call hdf5_read_attr(gr_halo_id, name, halopartNB)

             If(varlist%pos) Then
                Allocate(pos(3,halopartNB))
                name = 'position_part'
                Call hdf5_read_data(gr_halo_id, name, 3, halopartNB, pos)
             End If

             If(varlist%vel) Then
                Allocate(vel(3,halopartNB))
                name = 'velocity_part'
                Call hdf5_read_data(gr_halo_id, name, 3, halopartNB, vel)
             End If

             If(varlist%id) Then
                Allocate(id(halopartNB))
                name = 'identity_part'
                Call hdf5_read_data(gr_halo_id, name, halopartNB, id)
             End If

             If(varlist%pot) Then
                Allocate(pot(halopartNB))
                name = 'potential_part'
                Call hdf5_read_data(gr_halo_id, name, halopartNB, pot)
             End If

             Call hdf5_close_group(gr_halo_id)
             Call myfunction()

             If(allocated(id))  Deallocate(id)
             If(allocated(pos)) Deallocate(pos)
             If(allocated(vel)) Deallocate(vel)
             If(allocated(pot)) Deallocate(pot) 
             ih = ih + 1
          End Do
       End If
       ! next file 
       ifile = ifile + 1
       Deallocate(haloID)

       Call hdf5_close_file(file_id)
    End Do


  End Subroutine analyzelist

  !=======================================================================
  !> Analyze each halo contained in the halo file halofilename 
  !! Apply the function generic_function to each halo individually
  Subroutine analyzefile(halofilename, myfunction, inputvarlist)

    Use modiocommons
    
    Character(len=400), intent(in) :: halofilename
    Procedure(generic_function) :: myfunction
    Type(ListOfVar), intent(in), optional :: inputvarlist

    Type(ListOfVar)     :: varlist
    !Type(Type_common_metadata) :: common_metadata
    Integer :: i, ih, ifile            ! loop indices
    Integer :: ff, lf
    
    Namelist / Size / filelistsize    ! namelist
    Namelist / List / filelist        ! namelist
    
    If(Present(inputvarlist)) Then
       varlist = inputvarlist
    End If

#ifdef WITHMPI
    Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
    Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
    
    If(procID == 0) Then
#endif

       ! We read the list of halo files that we want to read and analyze: their number and name
       Open(Unit=10,file=trim(halofilename),status='old')
       Read(10,nml=Size)
       Allocate(filelist(filelistsize))
       Read(10,nml=List)
       Close(10)
       
       
       ff = 1
       lf = filelistsize
       myfilelistsize = filelistsize

#ifdef WITHMPI
    End If
    
    Call Mpi_Bcast(filelistsize, 1, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
    If(procID /= 0) Then
       Allocate(filelist(filelistsize))
    End If

    myfilelistsize = filelistsize / procNB
    ff = procID*myfilelistsize + 1
    If( procID < mod(filelistsize, procNB) ) Then
       myfilelistsize = myfilelistsize + 1
       ff = ff + procID
    Else
       ff = ff + mod(filelistsize,procNB)
    End If

    lf = ff + myfilelistsize - 1
    !print*,'procID,ff,lf:',procID,ff,lf

    Call Mpi_Bcast(filelist, len(filelist(1))*filelistsize, Mpi_Char, 0, Mpi_Comm_World, mpierr)
#endif

    ! Read param file and the Ramses Info Metadata.
    ! Also read Ramses cosmology file to compute hubble at z.
    ! For each processor, all this needs to be done ONLY ONCE.
    Do ifile = ff, ff

       ! open, read, close the param file
       open(17,file=trim(infile),form='formatted',status='old')
       read(17,*) box_len_Mpc  !bug of too restrictive format corrected by ry 20/02/2017
       read(17,*) simu_info%delta_vir !bug of too restrictive format corrected by ry 20/02/2017
       read(17,'(A)')    file_cosmo
       close(17)

       Call hdf5_open_file(filelist(ifile), file_id)
       Call h5read_meta_common(file_id, common_metadata)
       Call h5read_meta_pfof_parameter(file_id, param_pfof)
       Call h5read_meta_info_ramses(file_id, simu_info) 
       Call hdf5_close_file(file_id)

       !Definition of the convergence factor
       conv=(simu_info%h0/100./simu_info%aexp)

       !Definition of the mass of a particle, and the mean matter density of the universe
       !Assume not a zoom simulation. Need d0 and not integer to adapt to >2048**3 RY
       mass_part=1.d0/((2.d0**simu_info%levelmin)**simu_info%ndim)
       Rho_m=simu_info%unit_d

       !Critical density at redshift z=0 (in g/cm^3)
       Rho_crit0 = Rho_crit0 * (simu_info%h0/100.)**2

       !Compute hubble expansion at z (in a^2 H/H0 units)
       !Read ramses input file for cosmology (for hubble(z))
       open(90,file=trim(file_cosmo),status='old')
       NN=0
       DO 5 i=1,NTOT
         READ(90,*,END=6) ATEMP,EHUBBLETEMPA2,DUMMY,DUMMY,RTEMP
         AEXP_TEMP(NN+1)     = ATEMP
         EHUBBLE_TEMP(NN+1)  = EHUBBLETEMPA2/ATEMP**2
         RDIST_TEMP(NN+1)    = abs(RTEMP*CLIGHT/100.)
         NN=NN+1
5      CONTINUE
6      CLOSE(90)
       !Reverse order
       DO i=1,NN
        AEXP(i)     = AEXP_TEMP(NN-i+1)
        EHUBBLE(i)  = EHUBBLE_TEMP(NN-i+1)
        RDIST(i)    = RDIST_TEMP(NN-i+1)
       ENDDO

       !Find neighboring times
       i=3 !!!!!!WARNING: ramses_input.dat is not accurate enough very close to a=1 so I start at i=3 instead of 2!!!
       do while(AEXP(i)>simu_info%aexp .and. i<NN)
          i=i+1
       enddo
       !Interpolate
       Hz  = EHUBBLE(i)*(simu_info%aexp-AEXP(i-1))/(AEXP(i)-AEXP(i-1))+ &
           & EHUBBLE(i-1)*(simu_info%aexp-AEXP(i))/(AEXP(i-1)-AEXP(i))

       convdens=Hz**2/(simu_info%omega_m/simu_info%aexp**3) !rhoc/rhom

       Hz  = Hz * simu_info%aexp**2
       
       !Set delta_matter and deltacritz
       delta_matter=simu_info%delta_vir
       delta_critz=delta_matter/convdens
       If(procID==0) then
          print*,'Hz=',hz
          print*,'convdens',convdens
          print*,'aexp',simu_info%aexp
          print*,'delta_matter',delta_matter
          print*,'delta_critz',delta_critz
       endif
    EndDo

    If(procID==0) then
     print*, trim(infile)
     print*, box_len_Mpc
     print*, simu_info%delta_vir
     print*, file_cosmo

     print*, simu_info%ncpu
     print*, simu_info%ndim
     print*, simu_info%levelmin
     print*, simu_info%levelmax
     print*, simu_info%ngridmax
     print*, simu_info%nstep_coarse
     print*, simu_info%boxlen
     print*, simu_info%time
     print*, simu_info%aexp
     print*, simu_info%h0
     print*, simu_info%omega_m
     print*, simu_info%omega_l
     print*, simu_info%omega_k
     print*, simu_info%omega_b
     print*, simu_info%unit_l
     print*, simu_info%unit_d
     print*, simu_info%unit_t

     print *,'rho_crit0=',Rho_crit0
     print *,'Hubble(z)=',Hz
    Endif

    nbhaloanalyzed = 0
    ! We must first read the attribute nhalo_file from each file
    Do ifile = ff, lf
       Call hdf5_open_file(filelist(ifile), file_id)
       groupname = 'metadata'
       Call hdf5_open_group(file_id, groupname, gr_id)
       name = 'nhalo_file'
       Call hdf5_read_attr(gr_id, name, haloNB)
       if(procID==0)Print *,'nhalo_file=',ifile,haloNB
       nbhaloanalyzed = nbhaloanalyzed + haloNB
       Call hdf5_close_group(gr_id)
       Call hdf5_close_file(file_id)
    End Do

    currenthalo = 0

    Do ifile = ff, lf
       Call hdf5_open_file(filelist(ifile), file_id)
       groupname = 'metadata'
       Call hdf5_open_group(file_id, groupname, gr_id)
       name = 'nhalo_file'
       Call hdf5_read_attr(gr_id, name, haloNB)
       Call hdf5_close_group(gr_id)

       ! Read Pfof halo data
       groupname = 'data'
       Call hdf5_open_group(file_id, groupname, gr_id)
       name = 'identity_halo'
       Allocate(haloID(haloNB))
       Call hdf5_read_data(gr_id, name, haloNB, haloID)

       if(procID==0)print*,'ifile,haloNB:', ifile, haloNB
       Do ih = 1, haloNB
          
          currenthaloID = haloID(ih)
          !print*,'currenthaloID, haloID(ih)',currenthaloID, haloID(ih)
          currenthalo = currenthalo + 1
          
          groupname = "halo_0000000000000000000"
          Write(groupname(6:24),'(I19.19)') currenthaloID
          Call hdf5_open_group(gr_id, groupname, gr_halo_id)
          name = 'npart_halo'
          Call hdf5_read_attr(gr_halo_id, name, halopartNB)
          !print*,'npart_halo:',trim(groupname),halopartNB
          
          If(varlist%pos) Then
             Allocate(pos(3,halopartNB))
             name = 'position_part'
             Call hdf5_read_data(gr_halo_id, name, 3, halopartNB, pos)
          End If
          
          If(varlist%vel) Then
             Allocate(vel(3,halopartNB))
             name = 'velocity_part'
             Call hdf5_read_data(gr_halo_id, name, 3, halopartNB, vel)
          End If
          
          If(varlist%id) Then
             Allocate(id(halopartNB))
             name = 'identity_part'
             Call hdf5_read_data(gr_halo_id, name, halopartNB, id)
          End If

          If(varlist%pot) Then
             Allocate(pot(halopartNB))
             name = 'potential_part'
             Call hdf5_read_data(gr_halo_id, name, halopartNB, pot)
          End If
          
          Call hdf5_close_group(gr_halo_id)
          Call myfunction()
          
          If(allocated(id))  Deallocate(id)
          If(allocated(pos)) Deallocate(pos)
          If(allocated(vel)) Deallocate(vel)
          If(allocated(pot)) Deallocate(pot)
       End Do
       Deallocate(haloID)
       Call hdf5_close_group(gr_id)
       Call hdf5_close_file(file_id)
    End Do
  End Subroutine analyzefile

End Module modanalyzeengines
