!> @file
!! This file contains routines to write HDF5 halo files created with pfof and pfof_cone.

!> This module contains routines to write HDF5 halo files created with pfof and pfof_cone.
!>
!> Authors: F. Roy

Module modwritehalo

Contains

  !=======================================================================
  !> This subroutine writes for each halo the position, the velocity and the id of each particle in this halo
  !! in a hdf5 file.
  !! One file is written per MPI process and several halos per file.
  !! There is one group per halo, the name of the group is the ID of the halo.
  Subroutine h5writehalopart(mpicomm, haloNB, halopartNB, haloMass, haloID, halopartPos, halopartVel, &
       halopartID,halopartFor, halopartPot)

    Use modhdf5
    Use modparameters
    Use modmpicom
    Use modvariables
    Use modxdmf
    Implicit none
    
    Integer(kind=4), intent(in) :: haloNB !< Number of halos
    Integer(kind=4), intent(in) :: halopartNB  !< Number of particles belonging to the halos
    Integer(kind=PRI), dimension(haloNB), intent(in) :: haloID !< Array containing the ID of the halos
    Integer(kind=PRI), dimension(halopartNB), intent(in), target :: halopartID !< Array containing the ID of the
    !< particles belonging to the halos
    Integer(kind=4), dimension(haloNB), intent(in) :: haloMass !< Array containing the mass of the halos
    Real   (kind=4), dimension(3,halopartNB), intent(in), target :: halopartPos, halopartVel !< Array containing the position
    !< and velocity of the particles belonging to the halos
    Real   (kind=4), dimension(halopartNB), intent(in), target, optional :: halopartPot !< Array containing the potential of the
    !< particles belonging to the halos, this argument is optional
    Real   (kind=4), dimension(3,halopartNB), intent(in), target, optional :: halopartFor !< Array containing the force on the particles belonging to the halos, this argument is optional
    Integer(kind=4), intent(in) :: mpicomm !< MPI communicator used to gather informations 

    Integer(kind=4) :: ih, hptr
    Character(len=400) :: filestrct
    Character(len=400) :: filebase
    Character(len=5)  :: pid_char
    Character(len=16) :: groupname
    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=20) :: adata

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Integer(hid_t) :: gr_halo_id                            ! Group identifier

    Integer(kind=4), dimension(procNB) :: haloNBtab, dspl
    Integer(kind=4), dimension(:), allocatable :: haloMasstab
    Integer(kind=PRI), dimension(:), allocatable :: haloIDtab
    Integer(kind=PRI), dimension(2) :: IDminmax
    Integer(kind=4) :: haloNBall, p

    
#ifdef DEBUG
    Print *,"Enter h5writehalopart on process ",procID
#endif    

    Write(pid_char(1:5),'(I5.5)') procID
    filestrct = trim(output_root)//'_halo_'//pid_char//'.h5'
    filebase = trim(output_root)//'_halo'

    Call Mpi_Gather(haloNB, 1, Mpi_Integer, haloNBtab, 1, Mpi_Integer, 0, mpicomm, mpierr)
    If(procID==0) Then
       haloNBall = haloNBtab(1)
       dspl(1) = 0
       Do p=2, procNB 
          dspl(p) = haloNBall
          haloNBall= haloNBall + haloNBtab(p)
       End Do
       Allocate(haloMasstab(haloNBall), haloIDtab(haloNBall))
    Else
       Allocate(haloMasstab(0), haloIDtab(0))
    End If
    Call Mpi_Gatherv(haloMass, haloNB, Mpi_Integer, haloMasstab, haloNBtab, dspl, Mpi_Integer, 0, mpicomm, mpierr)
    Call Mpi_Gatherv(haloID, haloNB, Mpi_PRI, haloIDtab, haloNBtab, dspl, Mpi_PRI, 0, mpicomm, mpierr)

!!$    If(procID==0) Then
!!$       Call writesinglehaloxdmf(procNB,filebase,haloNBtab,haloMasstab,haloIDtab)
!!$    End If

    ! create the hdf5 file
    Call hdf5_create_file(filestrct, file_id, origin)
    
    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)
    
    aname = 'nfile'
    Call hdf5_write_attr(gr_root_id, aname, procNB)
    
    ! write the number of haloes as an attribute
    aname = "haloNB"
    Call hdf5_write_attr(gr_root_id, aname, haloNB)
    
    ! Write type as attribute    
    aname = 'Type'
    adata = 'halo format'
    Call hdf5_write_attr(gr_root_id, aname, adata)
    
    If(haloNB/=0) Then
       !! write the halo ID as data and not attribute: it seems that we cannot write integer(kind=8) attribute 
       aname = 'haloID'
       Call hdf5_write_data(gr_root_id, aname, haloNB, haloID)
       
       aname = 'haloIDminmax'
       IDminmax(1) = haloID(1)
       IDminmax(2) = haloID(haloNB)
       Call hdf5_write_data(gr_root_id, aname, 2, IDminmax)
       
       ! pointer to the current halo
       hptr = 1
       
       Do ih = 1, haloNB
          ! create a group for each halo
          groupname = "halo_00000000000"
          Write(groupname(6:16),'(I11.11)') haloID(ih)
          Call hdf5_create_group(gr_root_id, groupname, gr_halo_id)
          ! create an attribute containing the number of particles in the halo
          aname = "halopartNB"
          Call hdf5_write_attr(gr_halo_id, aname, haloMass(ih))
          
          dsetname="pos"
          Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartPos(:,hptr:hptr+haloMass(ih)-1)) 
          
          dsetname="vel"
          Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartVel(:,hptr:hptr+haloMass(ih)-1)) 
          
          dsetname = "ID"
          Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), halopartID(hptr:hptr+haloMass(ih)-1)) 
          
          If(present(halopartPot)) Then
             dsetname = 'pot'
             Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), halopartPot(hptr:hptr+haloMass(ih)-1)) 
          End If
          
          If(present(halopartFor)) Then
             dsetname = 'for'
             Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartFor(:,hptr:hptr+haloMass(ih)-1)) 
          End If
          
          ! Close the halo group.
          Call hdf5_close_group(gr_halo_id)
          
          ! move the pointer to the next halo
          hptr = hptr + haloMass(ih)
       End Do
       
       ! Close the root group.
       Call hdf5_close_group(gr_root_id)
       
       Call hdf5_close_file(file_id)
       
    End If

#ifdef DEBUG
    Print *,"Exit h5writehalopart on process ",procID
#endif

  End Subroutine h5writehalopart



  !=======================================================================
  !> This subroutine writes, for each halo, its mass (number of particles), 
  !! the position and the velocity of its center of mass and its ID
  !! in only one hdf5 file using parallel HDF5.
  Subroutine mpih5writehalomass(mpicomm, haloNB_all, haloNB,nh, haloMass, halocomPos, halocomVel, &
       haloID, haloRadius, haloSubHaloNB)

    Use modparameters
    Use modvariables
    Use modmpicom
    Use modhdf5

    Implicit none
    
    Integer(kind=4), intent(in) :: haloNB_all !< Total number of halos (sum over all processes)
    Integer(kind=4), intent(in) :: haloNB !< Number of halos
    Integer(kind=4), intent(in) :: nh 
!!$    Integer(kind=4), dimension(haloNB), intent(in), target :: haloMass !< Mass of the halos
!!$    Real(kind=8), dimension(3,haloNB), intent(in), target :: halocomPos !< Position of the center of mass of the halos
!!$    Real(kind=8), dimension(3,haloNB), intent(in), target :: halocomVel !< Velocity of the center of mass of the halos
!!$    Integer(kind=PRI), dimension(haloNB), intent(in), target :: haloID !< ID of the halos
!!$    Real(kind=8), dimension(haloNB), intent(in), target :: haloRadius !< Radius of the halos
!!$    Integer(kind=4), dimension(haloNB), intent(in), target :: haloSubHaloNB !< Number of subhalos in each halo
    Integer(kind=4), dimension(nh), intent(in), target :: haloMass !< Mass of the halos
    Real(kind=8), dimension(3,nh), intent(in), target :: halocomPos !< Position of the center of mass of the halos
    Real(kind=8), dimension(3,nh), intent(in), target :: halocomVel !< Velocity of the center of mass of the halos
    Integer(kind=PRI), dimension(nh), intent(in), target :: haloID !< ID of the halos
    Real(kind=8), dimension(nh), intent(in), target :: haloRadius !< Radius of the halos
    Integer(kind=4), dimension(nh), intent(in), target :: haloSubHaloNB !< Number of subhalos in each halo

    Integer(kind=4), intent(in) :: mpicomm !< MPI communicator used to create and write the file

    Character(len=400) :: filename
    Character(len=16) :: aname                           ! Attribute name
    Character(len=16) :: dsetname                        ! Dataset name
    Character(len=20) :: adata
    Integer(hid_t) :: file_id                            ! File identifier
    Integer(kind=4) :: IntSubHalo
    Integer(kind=4) :: IntUnbinding
    Integer(kind=4) :: begh, endh
    Logical(kind=4) :: empty

#ifdef DEBUG
    Print *,"Enter mpih5writehalomass on process ",procID
#endif    
    
    filename = trim(output_root)//'_halomass.h5'

    ! Create h5 parallel file
    Call hdf5_create_mpi_file(filename, mpicomm, file_id, origin)

    ! Write number of halos as attribute
    aname = 'haloNB'
    Call hdf5_write_attr(file_id, aname, haloNB_all)

    ! Write type as attribute    
    aname = 'Type'
    adata = 'halomass format'
    Call hdf5_write_attr(file_id, aname, adata)

    ! Write percolation parameter b as attribute
    aname = 'perco'
    Call hdf5_write_attr(file_id, aname, perco)

    ! Write minimum halo mass Mmin as attribute
    aname = 'Mmin'
    Call hdf5_write_attr(file_id, aname, Mmin)

    ! Write doUnbinding as attribute (not implemented yet) 
    aname = 'doUnbinding'
    If(doUnbinding) Then
       IntUnbinding=1
    Else
       IntUnbinding=0
    End If
    Call hdf5_write_attr(file_id, aname, IntUnbinding)

    ! Write doSubHalo as attribute (not implemented yet) 
    aname = 'doSubHalo'
    If(doSubHalo) Then
       IntSubHalo=1
    Else
       IntSubHalo=0
    End If
    Call hdf5_write_attr(file_id, aname, IntSubHalo)

    ! Write the simulation name as attribute
    aname = 'SimulationName'
    Call hdf5_write_attr(file_id, aname, simulationName)
    
    ! Write the output number as attribute (for standard pfof only)
    aname = 'outputNB'
    Call hdf5_write_attr(file_id, aname, outputNB)

    If(haloNB == 0) Then
!       Print *,'HALONB=0!!!!!',size(halocomPos,2)
       
!       nh = 1
       begh = 1
       endh = 1
       empty = .true.
    Else
       begh = 1
       endh = haloNB
       empty = .false.
    End If

    dsetname = 'halocomPos'
    Call hdf5_write_mpi_data(file_id, dsetname, 3, nh, halocomPos(:,begh:endh), mpicomm, empty)

    dsetname = 'halocomVel'
    Call hdf5_write_mpi_data(file_id, dsetname, 3, nh, halocomVel(:,begh:endh), mpicomm, empty)

    dsetname = 'haloID'
    Call hdf5_write_mpi_data(file_id, dsetname, nh, haloID(begh:endh), mpicomm, empty)

    dsetname = 'haloMass'
    Call hdf5_write_mpi_data(file_id, dsetname, nh, haloMass(begh:endh), mpicomm, empty)

    dsetname = 'haloRadius'
    Call hdf5_write_mpi_data(file_id, dsetname, nh, haloRadius(begh:endh), mpicomm, empty)

    If(doSubHalo) Then
       dsetname = 'haloSubHaloNB'
       Call hdf5_write_mpi_data(file_id, dsetname, nh, haloSubHaloNB(begh:endh), mpicomm, empty)
    End If

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

#ifdef DEBUG
    Print *,"Exit mpih5writehalomass on process ",procID
#endif    
    
  End Subroutine mpih5writehalomass


End Module modwritehalo
