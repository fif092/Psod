!==============================================================================
! Project: pFoF
! File: common/src/modreadinfo.f90
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
!! This file contains subroutines to read info files written by RAMSES for snapshot and light cone.

!> This module contains subroutines to read info files written by RAMSES for snapshot and light cone.
!>
!> Authors: F. Roy

Module modreadinfo

  Use modconstant

  Implicit none

  Private

  Public :: readinforamses, &
       readinfoconepart, &
       readinfoconegrav

Contains

  !=======================================================================
  !> Read RAMSES snapshot info file
  Subroutine readinforamses(filename, inforamses, ierr, errormessage)
    
    Implicit none

    Character(len=400), intent(in) :: filename        !< info file name
    Type(Type_info_ramses), intent(out) :: inforamses !< info structure
    Integer(kind=4), intent(out) :: ierr              !< error code
    Character(len=500), intent(out) :: errormessage   !< error message

    ! Local variable
    Character(len=13) :: dumchar
   

    Open(Unit=12, file=filename, status='old', iostat=ierr)
    If( ierr > 0 ) Then
       errormessage='Error opening ramses info file '//trim(filename)
       Return
    End If
    
    Read(12,'(A13,I11)', iostat=ierr) dumchar, inforamses%ncpu
    Read(12,'(A13,I11)', iostat=ierr) dumchar, inforamses%ndim
    Read(12,'(A13,I11)', iostat=ierr) dumchar, inforamses%levelmin
    Read(12,'(A13,I11)', iostat=ierr) dumchar, inforamses%levelmax
    Read(12,'(A13,I11)', iostat=ierr) dumchar, inforamses%ngridmax
    Read(12,'(A13,I11)', iostat=ierr) dumchar, inforamses%nstep_coarse
    Read(12,*, iostat=ierr) 
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%boxlen
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%time
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%aexp
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%h0
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%omega_m
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%omega_l
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%omega_k
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%omega_b
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%unit_l
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%unit_d
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, inforamses%unit_t
    Close(12)

    If( ierr > 0 ) Then
       errormessage='Error reading ramses info file '//trim(filename)
       Return
    End If    
    
  End Subroutine readinforamses


  !=======================================================================
  !> Read RAMSES particles light cone info file
  Subroutine readinfoconepart(filename, infocone, ierr, errormessage)
    
    Implicit none

    Character(len=400), intent(in) :: filename         !< info file name
    Type(Type_info_cone_part), intent(out) :: infocone !< info structure
    Integer(kind=4), intent(out) :: ierr               !< error code
    Character(len=500), intent(out) :: errormessage    !< error message

    ! Local variable
    Character(len=13) :: dumchar

    Open(Unit=12, file=filename, status='old', iostat=ierr)
    If(ierr > 0) Then
       errormessage = 'Error opening conepart info file '//trim(filename)
       Return
    End If
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%ncpu
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%nstride
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%nstep_coarse
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aexp
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_x
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_y
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_z
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_rds
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%cone_id
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%cone_zlim
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%amax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%amin
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zmax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zmin
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dmax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dmin
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dtol
    Read(12,'(A13,I20)', iostat=ierr) dumchar, infocone%nglobalfile
    Read(12,'(A13,I20)', iostat=ierr) dumchar, infocone%npart
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%isfullsky
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%thetay
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%thetaz
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%theta
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%phi
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aendconem2
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aendconem1
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aendcone
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aexpold
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aexp
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zendconem2
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zendconem1
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zendcone
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zexpold
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zexp
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dendconem2
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dendconem1
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dendcone
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dexpold
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dexp
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%future              
    Close(12)

    If(ierr > 0) Then
       errormessage = 'Error reading in conepart info file '//trim(filename)
       Return
    End If

  End Subroutine readinfoconepart

  !=======================================================================
  !> Read RAMSES cells light cone info file
  Subroutine readinfoconegrav(filename, infocone, ierr, errormessage)
    
    Implicit none

    Character(len=400), intent(in) :: filename         !< info file name
    Type(Type_info_cone_grav), intent(out) :: infocone !< info structure
    Integer(kind=4), intent(out) :: ierr               !< error code
    Character(len=500), intent(out) :: errormessage    !< error message

    ! Local variable
    Character(len=13) :: dumchar

    Open(Unit=12, file=filename, status='old', iostat=ierr)
    If(ierr > 0) Then
       errormessage = 'Error opening conegrav info file '//trim(filename)
       Return
    End If
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%ncpu
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%nstride
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%nstep_coarse
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aexp
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%nlevel
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%levelmin
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%levelmax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_x
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_y
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_z
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%observer_rds
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%cone_id
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%cone_zlim
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%amax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%amin
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zmax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zmin
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dmax
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dmin
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dtol
    Read(12,'(A13,I20)', iostat=ierr) dumchar, infocone%nglobalfile
    Read(12,'(A13,I20)', iostat=ierr) dumchar, infocone%nglobalcell
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%isfullsky
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%thetay
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%thetaz
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%theta
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%phi
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aendconem2
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aendconem1
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%aendcone
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zendconem2
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zendconem1
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%zendcone
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dendconem2
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dendconem1
    Read(12,'(A13,E24.15)', iostat=ierr) dumchar, infocone%dendcone
    Read(12,'(A13,I11)', iostat=ierr) dumchar, infocone%future
    Close(12)

    If(ierr > 0) Then
       errormessage = 'Error reading in conegrav info file '//trim(filename)
       Return
    End If

  End Subroutine readinfoconegrav

End Module modreadinfo
