!==============================================================================
! Project: pFoF
! File: common/src/modwritehalo.f90
! Copyright Fabrice Roy (2015)
! Fabrice.Roy@obspm.fr
! 
! This software is a computer program whose purpose is to detect dark matter
! haloes in cosmological simulation with a parallel Friend of Friend algorithm.
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


!> @file
!! This file contains routines to write HDF5 halo files created with pfof and pfof_cone.

!> This module contains routines to write HDF5 halo files created with pfof and pfof_cone.
!>
!> Authors: F. Roy

Module modwritehalo

  Use modiocommons
  Use modconstant

  Implicit none

  Private

  Public :: h5writehalopart, &
       mpih5writehalopart, &
       mpih5writehalomass

Contains

  !=======================================================================
  !> This subroutine writes for each halo the position, the velocity and the id of each particle in this halo
  !! in a hdf5 file.
  !! One file is written per MPI process and several halos per file.
  !! There is one group per halo, the name of the group is the ID of the halo.
  Subroutine h5writehalopart(info_proc, param_pfof, haloNB, halopartNB, haloMass, haloID, &
       halopartPos, halopartVel, halopartID, halopartFor, halopartPot, halopartRamsesID, &
       inforamses, infocone)

    Use modhdf5
    Use mpi
    Use modmpicommons !, only : procID, procNB

    Implicit none
    
    Type(Type_info_process), intent(in) :: info_proc !< MPI communicator used to gather informations 
    Class(Type_parameter_halofinder), intent(in) :: param_pfof
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
    Integer(kind=PRI), dimension(halopartNB), intent(in), target, optional :: halopartRamsesID !< Array containing the "Ramses" ID of the particles belonging to the halos detected in lightcones

    Type(Type_info_ramses),intent(in) :: inforamses
    Class(Type_info_cone), intent(in), optional :: infocone

    Integer(kind=4) :: ih, hptr
    Character(len=400) :: filestrct
    Character(len=5)  :: pid_char
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: halofinder
    
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    Integer(hid_t) :: gr_halo_id                            ! Group identifier

    Integer(kind=4), dimension(procNB) :: haloNBtab, dspl
    Integer(kind=4), dimension(:), allocatable :: haloMasstab
    Integer(kind=PRI), dimension(:), allocatable :: haloIDtab
    Integer(kind=PRI), dimension(2) :: IDminmax
    Integer(kind=4) :: haloNBall, p
    Integer(kind=4) :: mpierr

    Integer(kind=8) :: npart8
    Character(len=H5STRLEN) :: codename
    Logical(kind=4) :: islast

#ifdef DEBUG
    Print *,'Enter h5writehalopart on process ', info_proc%global_comm%pid
#endif    

    Write(pid_char(1:5),'(I5.5)')  info_proc%global_comm%pid

    Select Type (param_pfof)
    Class is (Type_parameter_pfof_snap)
       filestrct = 'pfof_halo_snap_part_data_'//&
            trim(param_pfof%simulation_name)//'_'//pid_char//'.h5'
       codename='pfof_snap'
       halofinder = 'pfof'
       islast = .false.
    Class is (Type_parameter_pfof_cone)
       filestrct = 'pfof_halo_cone_part_data_'//&
            trim(param_pfof%simulation_name)//'_'//pid_char//'.h5'
       codename='pfof_cone'
       halofinder = 'pfof'
       islast = .true.
    Class is (Type_parameter_psod_snap)
       filestrct = 'psod_halo_snap_part_data_'//&
            trim(param_pfof%simulation_name)//'_'//pid_char//'.h5'
       codename='psod_snap'
       halofinder = 'psod'
       islast = .false.
    End Select


    Call Mpi_Gather(haloNB, 1, Mpi_Integer, haloNBtab, 1, Mpi_Integer, 0, &
         info_proc%global_comm%name, mpierr)
    If(info_proc%global_comm%pid == 0) Then
       haloNBall = haloNBtab(1)
       dspl(1) = 0
       Do p=2, info_proc%global_comm%size
          dspl(p) = haloNBall
          haloNBall= haloNBall + haloNBtab(p)
       End Do
       Allocate(haloMasstab(haloNBall), haloIDtab(haloNBall))
    Else
       Allocate(haloMasstab(0), haloIDtab(0))
    End If
    Call Mpi_Gatherv(haloMass, haloNB, Mpi_Integer, haloMasstab, haloNBtab, dspl,&
         Mpi_Integer, 0, info_proc%global_comm%name, mpierr)
    Call Mpi_Gatherv(haloID, haloNB, Mpi_PRI, haloIDtab, haloNBtab, dspl, &
         Mpi_PRI, 0, info_proc%global_comm%name, mpierr)

    ! create the hdf5 file
    Call hdf5_create_file(filestrct, file_id)

    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8, info_proc%global_comm%pid)
    Call h5write_meta_halofinder_parameter(file_id, param_pfof)
    Call h5write_meta_info_ramses(file_id, inforamses, islast)

    If(present(infocone)) Then
       Call h5write_meta_info_cone(file_id, infocone,islast)
    End If


    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    aname = 'halo_finder'
    adata = halofinder
    Call hdf5_write_attr(gr_id, aname, adata)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'halo'
    Call hdf5_write_attr(gr_id, aname, adata)
    
    aname = 'nfile'
    Call hdf5_write_attr(gr_id, aname, info_proc%global_comm%size)

    ! write the number of haloes as an attribute
    aname = 'nhalo_file'
    Call hdf5_write_attr(gr_id, aname, haloNB)

    ! write the number of particles written in the file
    dsetname = 'npart_file'
    npart8 = int(halopartNB,kind=8)
    Call hdf5_write_data(gr_id, dsetname, npart8)

    Call hdf5_close_group(gr_id)
    
    groupname = 'data'
    Call hdf5_create_group(file_id, groupname, gr_id)

    If(haloNB/=0) Then
       !! write the halo ID as data and not attribute: it seems that we cannot write integer(kind=8) attribute 
       aname = 'identity_halo'
       Call hdf5_write_data(gr_id, aname, haloNB, haloID)
       
       aname = 'identity_halo_minmax'
       IDminmax(1) = haloID(1)
       IDminmax(2) = haloID(haloNB)
       Call hdf5_write_data(gr_id, aname, 2, IDminmax)
       
       ! pointer to the current halo
       hptr = 1
       
       Do ih = 1, haloNB
          ! create a group for each halo
          groupname = 'halo_0000000000000000000'
          Write(groupname(6:24),'(I19.19)') haloID(ih)
          Call hdf5_create_group(gr_id, groupname, gr_halo_id)
          ! create an attribute containing the number of particles in the halo
          aname = 'npart_halo'
          Call hdf5_write_attr(gr_halo_id, aname, haloMass(ih))
          
          dsetname='position_part'
          Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), &
               halopartPos(:,hptr:hptr+haloMass(ih)-1)) 
          
          dsetname='velocity_part'
          Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), &
               halopartVel(:,hptr:hptr+haloMass(ih)-1)) 
          
          dsetname = 'identity_part'
          Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), &
               halopartID(hptr:hptr+haloMass(ih)-1)) 
          
          If(present(halopartPot)) Then
             dsetname = 'potential_part'
             Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), &
                  halopartPot(hptr:hptr+haloMass(ih)-1)) 
          End If
          
          If(present(halopartFor)) Then
             dsetname = 'gravitational_field_part'
             Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), &
                  halopartFor(:,hptr:hptr+haloMass(ih)-1)) 
          End If
          
          If(present(halopartRamsesID)) Then
             dsetname = 'ramses_identity_part'
             Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), &
                  halopartRamsesID(hptr:hptr+haloMass(ih)-1))
          End If

          ! Close the halo group.
          Call hdf5_close_group(gr_halo_id)
          
          ! move the pointer to the next halo
          hptr = hptr + haloMass(ih)
       End Do
       
       ! Close the root group.
       Call hdf5_close_group(gr_id)
       
       Call hdf5_close_file(file_id)
       
    End If

#ifdef DEBUG
    Print *,'Exit h5writehalopart on process ',info_proc%global_comm%pid
#endif

  End Subroutine h5writehalopart



  !=======================================================================
  !> This subroutine writes for each halo the position, the velocity and the id of each
  !! particle in this halo in a hdf5 file.
  !! One file is written per MPI communicator write_comm and several halos per file, using parallel HDF5.
  !! There is one group per halo, the name of the group is the ID of the halo.
  Subroutine mpih5writehalopart(info_proc, param_pfof, haloNB, halopartNB, haloMass, haloID, &
       halopartPos, halopartVel, halopartID, halopartFor, halopartPot, halopartRamsesID, &
       inforamses, infocone)

    Use modhdf5
    Use mpi
    Use modmpicommons !, only : procID, procNB

    Implicit none
    
    Type(Type_info_process), intent(in) :: info_proc !< MPI communicator used to gather informations 
    Class(Type_parameter_halofinder), intent(in) :: param_pfof
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
    Integer(kind=PRI), dimension(halopartNB), intent(in), target, optional :: halopartRamsesID !< Array containing the "Ramses" ID of the particles belonging to the halos detected in lightcones

    Type(Type_info_ramses),intent(in) :: inforamses
    Class(Type_info_cone), intent(in), optional :: infocone

    Integer(kind=4) :: ih
    Character(len=400) :: filestrct
    Character(len=5)  :: pid_char
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: halofinder
    
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_halo_id                            ! Group identifier
    Integer(hid_t) :: gr_data_id                            ! Group identifier
    Integer(hid_t) :: gr_meta_id                            ! Group identifier
    
    Integer(kind=4), dimension(:), allocatable :: haloNBtab, dspl
    Integer(kind=4), dimension(:), allocatable :: haloMasstab
    Integer(kind=PRI), dimension(:), allocatable :: haloIDtab
    Integer(kind=PRI), dimension(2) :: IDminmax
    Integer(kind=4) :: haloNBall, p

    Integer(kind=8) :: npart8
    Character(len=H5STRLEN) :: codename
    Logical(kind=4) :: islast

    Integer(kind=4) :: fp, lp, hptr, pptr, partnb
    Logical(kind=4) :: empty

    Integer(kind=4) :: mpierr
    Integer(kind=4) :: nfile
    Integer(kind=4), dimension(procNB) :: colortab

#ifdef DEBUG
    Print *,'Enter mpih5writehalopart on process ',procID
#endif    

    ! pid_char contains the 'index' of the file in which the process will write its data
    Write(pid_char(1:5),'(I5.5)') info_proc%write_comm%color

    Allocate(haloNBtab(info_proc%write_comm%size), dspl(info_proc%write_comm%size))
    
    Select Type (param_pfof)
    Class is (Type_parameter_pfof_snap)
       filestrct = 'pfof_halo_snap_part_data_'//&
            trim(param_pfof%simulation_name)//'_'//pid_char//'.h5'
       codename='pfof_snap'
       halofinder = 'pfof'
       islast = .false.
    Class is (Type_parameter_pfof_cone)
       filestrct = 'pfof_halo_cone_part_data_'//&
            trim(param_pfof%simulation_name)//'_'//pid_char//'.h5'
       codename='pfof_cone'
       halofinder = 'pfof'
       islast = .true.
    Class is (Type_parameter_psod_snap)
       filestrct = 'psod_halo_snap_part_data_'//&
            trim(param_pfof%simulation_name)//'_'//pid_char//'.h5'
       codename='psod_snap'
       halofinder = 'psod'
       islast = .false.
    End Select

    !! essai...
    Call Mpi_Allgather(info_proc%write_comm%color, 1, Mpi_Integer, colortab, 1, Mpi_Integer, &
         info_proc%global_comm%name, mpierr)
    nfile = maxval(colortab) - minval(colortab) + 1
    
    ! Gather halo nb from every process on process 0
    Call Mpi_Gather(haloNB, 1, Mpi_Integer, haloNBtab, 1, Mpi_Integer, 0, &
         info_proc%write_comm%name, mpierr)

   ! Create arrays containing halo mass and halo ID with length haloNBall 
    If(info_proc%write_comm%pid==0) Then
       haloNBall = haloNBtab(1)
       dspl(1) = 0
       Do p=2, info_proc%write_comm%size
          dspl(p) = haloNBall
          haloNBall= haloNBall + haloNBtab(p)
       End Do
    End If

    Call Mpi_Bcast(haloNBall, 1, Mpi_Integer, 0, info_proc%write_comm%name, mpierr) 
    Allocate(haloMasstab(haloNBall), haloIDtab(haloNBall))

    ! Process 0 gathers mass and ID of each halo
    Call Mpi_Gatherv(haloMass, haloNB, Mpi_Integer, haloMasstab, haloNBtab, dspl, &
         Mpi_Integer, 0, info_proc%write_comm%name, mpierr)

    Call Mpi_Gatherv(haloID, haloNB, Mpi_PRI, haloIDtab, haloNBtab, dspl, &
         Mpi_PRI, 0, info_proc%write_comm%name, mpierr)
    
    Call Mpi_Bcast(haloIDtab, haloNBall, MPI_PRI, 0, info_proc%write_comm%name, mpierr)

    Call Mpi_Bcast(haloMasstab, haloNBall, Mpi_Integer, 0, info_proc%write_comm%name, mpierr)
    
    ! create the hdf5 file
    Call hdf5_create_mpi_file(filestrct, info_proc%write_comm%name, file_id)

    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8, procID)
    Call h5write_meta_halofinder_parameter(file_id, param_pfof)
    Call h5write_meta_info_ramses(file_id, inforamses, islast)

    If(present(infocone)) Then
       Call h5write_meta_info_cone(file_id, infocone,islast)
    End If

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_meta_id)

    aname = 'halo_finder'
    adata = halofinder
    Call hdf5_write_attr(gr_meta_id, aname, adata)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'halo'
    Call hdf5_write_attr(gr_meta_id, aname, adata)
    
    aname = 'nfile'
    Call hdf5_write_attr(gr_meta_id, aname, nfile )

    ! write the number of haloes as an attribute
    aname = 'nhalo_file'
    Call hdf5_write_attr(gr_meta_id, aname, haloNBall)

    ! write the number of particles written in the file
    dsetname = 'npart_file'
    npart8 = int(sum(haloMasstab) ,kind=8)
    Call hdf5_write_data(gr_meta_id, dsetname, npart8)

    Call hdf5_close_group(gr_meta_id)

    groupname = 'data'
    Call hdf5_create_group(file_id, groupname, gr_data_id)
    
    If(haloNBall/=0) Then
          
       If(haloNB==0) Then
          empty=.true.
       Else
          empty=.false.
       End If

       !! write the halo ID as data and not attribute:
       !! it seems that we cannot write integer(kind=8) attribute
       !! parallel write: each process write its haloID array
       aname = 'identity_halo'
       Call hdf5_write_mpi_data(gr_data_id, aname, haloNB, haloID,&
            info_proc%write_comm%name,empty) 

       If(info_proc%write_comm%pid == 0) Then
          IDminmax(1) = haloIDtab(1)
          IDminmax(2) = haloIDtab(haloNBall)
          empty = .false.
       Else
          IDminmax(1) = haloID(1)
          IDminmax(2) = haloID(haloNB)
          empty = .true.
       End If
       aname = 'identity_halo_minmax'
       Call hdf5_write_mpi_data(gr_data_id, aname, 2, IDminmax, &
            info_proc%write_comm%name, empty)

       ! pointer to the current halo
       hptr = 1
       ! pointer to the current particle
       pptr = 1
       
       Do ih = 1, haloNBall
          ! create a group for each halo
          groupname = 'halo_0000000000000000000'
          Write(groupname(6:24),'(I19.19)') haloIDtab(ih)
          Call hdf5_create_group(gr_data_id, groupname, gr_halo_id)

          If( haloNB == 0 .or. hptr > haloNB) Then
             empty = .true.
             fp = 1
             lp = 1
          Else
             If( haloID(hptr) == haloIDtab(ih) ) Then
                empty = .false.
                partnb = haloMass(hptr)
                fp = pptr
                lp = pptr + partnb - 1
                hptr = hptr + 1
                pptr = pptr + partnb
             Else
                empty = .true.
                partnb = 1
                fp = 1
                lp = 1
             End If
          End If
          
          ! create an attribute containing the number of particles in the halo
          aname = 'npart_halo'
          Call hdf5_write_attr(gr_halo_id, aname, haloMasstab(ih))
          
          dsetname='position_part'
          Call hdf5_write_mpi_data(gr_halo_id, dsetname, 3, partnb, &
               halopartPos(:,fp:lp), info_proc%write_comm%name, empty) 

          
          dsetname='velocity_part'
          Call hdf5_write_mpi_data(gr_halo_id, dsetname, 3, partnb, &
               halopartVel(:,fp:lp), info_proc%write_comm%name, empty) 

          dsetname = 'identity_part'
          Call hdf5_write_mpi_data(gr_halo_id, dsetname, partnb, &
               halopartID(fp:lp), info_proc%write_comm%name, empty) 
          
          If(present(halopartPot)) Then
             dsetname = 'potential_part'
             Call hdf5_write_mpi_data(gr_halo_id, dsetname, partnb, &
                  halopartPot(fp:lp), info_proc%write_comm%name, empty) 
          End If
          
          If(present(halopartFor)) Then
             dsetname = 'gravitational_field_part'
             Call hdf5_write_mpi_data(gr_halo_id, dsetname, 3, partnb, &
                  halopartFor(:,fp:lp), info_proc%write_comm%name, empty) 
          End If
          
          If(present(halopartRamsesID)) Then
             dsetname = 'ramses_identity_part'
             Call hdf5_write_mpi_data(gr_halo_id, dsetname, partnb, &
                  halopartRamsesID(fp:lp), info_proc%write_comm%name, empty)
          End If
          
          ! Close the halo group.
          Call hdf5_close_group(gr_halo_id)
                    
       End Do

    End If


    Call hdf5_close_group(gr_data_id)
    Call hdf5_close_mpi_file(file_id)

    Deallocate(haloNBtab, dspl, haloMasstab, haloIDtab)
    
#ifdef DEBUG
    Print *,'Exit mpih5writehalopart on process ',procID
#endif

  End Subroutine mpih5writehalopart



  !=======================================================================
  !> This subroutine writes, for each halo, its mass (number of particles), 
  !! the position and the velocity of its center of mass and its ID
  !! in only one hdf5 file using parallel HDF5.
  Subroutine mpih5writehalomass(mpicomm, param_pfof, haloNB_all, haloNB, nh, &
       haloMass, halocomPos, halocomVel, haloID, haloRadius, haloSubHaloNB, &
       inforamses, infocone)

    Use modmpicommons, only : procID, EmergencyStop
    Use modhdf5

    Implicit none

    Integer(kind=4), intent(in) :: haloNB_all !< Total number of halos (sum over all processes)
    Integer(kind=4), intent(in) :: haloNB !< Number of halos
    Integer(kind=4), intent(in) :: nh 
    Integer(kind=4), dimension(nh), intent(in), target :: haloMass !< Mass of the halos
    Real(kind=8), dimension(3,nh), intent(in), target :: halocomPos !< Position of the center of mass of the halos
    Real(kind=8), dimension(3,nh), intent(in), target :: halocomVel !< Velocity of the center of mass of the halos
    Integer(kind=PRI), dimension(nh), intent(in), target :: haloID !< ID of the halos
    Real(kind=8), dimension(nh), intent(in), target :: haloRadius !< Radius of the halos
    Integer(kind=4), dimension(nh), intent(in), target :: haloSubHaloNB !< Number of subhalos in each halo
    Integer(kind=4), intent(in) :: mpicomm !< MPI communicator used to create and write the file
    Type(Type_info_ramses),intent(in) :: inforamses
    Class(Type_info_cone), intent(in), optional :: infocone
    Class(Type_parameter_pfof), intent(in) :: param_pfof

    Character(len=400) :: filename
    Character(len=H5STRLEN) :: aname                           ! Attribute name
    Character(len=H5STRLEN) :: dsetname                        ! Dataset name
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: codename
    Character(len=H5STRLEN) :: halofinder
    Character(len=400) :: errormessage

    Integer(hid_t) :: file_id                            ! File identifier
    Integer(hid_t) :: gr_id
    Integer(kind=4) :: begh, endh
    Logical(kind=4) :: empty
    Integer(kind=8) :: npart8
    Logical(kind=4) :: islast

#ifdef DEBUG
    Print *,'Enter mpih5writehalomass on process ',procID
#endif    
    
    Select Type(param_pfof)
    Class Is(Type_parameter_pfof_snap)
       filename = 'pfof_halo_snap_part_hfprop_'//trim(param_pfof%simulation_name)//'.h5'
       codename = 'pfof_snap'
       halofinder = 'pfof'
       islast = .false.
    Class Is(Type_parameter_pfof_cone)
       filename = 'pfof_halo_cone_part_hfprop_'//trim(param_pfof%simulation_name)//'.h5'
       codename = 'pfof_cone'
       halofinder = 'pfof'
       islast = .true.
    End Select

#ifdef DEBUG
    Print *,'Enter mpih5writehalomass: filename is ',trim(filename),' on process ',procID
#endif    
    
    ! Create h5 parallel file
    Call hdf5_create_mpi_file(filename, mpicomm, file_id)

    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8, procID)
    Call h5write_meta_pfof_parameter(file_id, param_pfof)
    Call h5write_meta_info_ramses(file_id, inforamses,islast)

    If(present(infocone)) Then
       Call h5write_meta_info_cone(file_id, infocone,islast)
    End If


    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    aname = 'halo_finder'
    adata = halofinder
    Call hdf5_write_attr(gr_id, aname, adata)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'halomass'
    Call hdf5_write_attr(gr_id, aname, adata)
    
    aname = 'center_type'
    adata = 'center_of_mass'
    Call hdf5_write_attr(gr_id, aname, adata)
    
    ! write the number of haloes as an attribute
    aname = 'nhalo_file'
    Call hdf5_write_attr(gr_id, aname, haloNB_all)

    Call hdf5_close_group(gr_id)

    If(haloNB == 0) Then
       begh = 1
       endh = 1
       empty = .true.
    Else
       begh = 1
       endh = haloNB
       empty = .false.
    End If

#ifdef DEBUG
    Print *,procID, 'Sizes of arrays in mpih5writehalomass:',nh,endh-begh+1
    If(nh/=endh-begh+1) Then
       errormessage = 'Error in mpih5writehalomass'
       Call EmergencyStop(errormessage,111)
    End If
#endif
    
    groupname = 'data'
    Call hdf5_create_group(file_id, groupname, gr_id)

    dsetname = 'position_halo'
    Call hdf5_write_mpi_data(gr_id, dsetname, 3, nh, halocomPos(:,begh:endh), mpicomm, empty)

    dsetname = 'velocity_halo'
    Call hdf5_write_mpi_data(gr_id, dsetname, 3, nh, halocomVel(:,begh:endh), mpicomm, empty)

    dsetname = 'identity_halo'
    Call hdf5_write_mpi_data(gr_id, dsetname, nh, haloID(begh:endh), mpicomm, empty)

    dsetname = 'npart_halo'
    Call hdf5_write_mpi_data(gr_id, dsetname, nh, haloMass(begh:endh), mpicomm, empty)

    dsetname = 'rmax_halo'
    Call hdf5_write_mpi_data(gr_id, dsetname, nh, haloRadius(begh:endh), mpicomm, empty)
    
    If(param_pfof%do_subhalo) Then
       dsetname = 'nsubhalo_halo'
       Call hdf5_write_mpi_data(gr_id, dsetname, nh, haloSubHaloNB(begh:endh), mpicomm, empty)
    End If

    Call hdf5_close_group(gr_id)

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

#ifdef DEBUG
    Print *,'Exit mpih5writehalomass on process ',procID
#endif    
    
  End Subroutine mpih5writehalomass


End Module modwritehalo
