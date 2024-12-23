Module modreadinfo
  
  Use modconstant
  Use modmpicommons

  Implicit none

Contains

  !=======================================================================
  Subroutine readinforamses(filename, procID, inforamses)
    
    Implicit none

    Character(len=400), intent(in) :: filename
    Integer(kind=4), intent(in) :: procID
    Type(Type_inforamses), intent(out) :: inforamses

    ! Local variable
    Integer(kind=4) :: ierr
    Character(len=13) :: dumchar

    If(procID==0) Then
       Open(Unit=12, file=filename, status='old', iostat=ierr)
       If(ierr/=0) Then
          Print *,'Error opening ramses info file ',filename
       End If

       Read(12,'(A13,I11)') dumchar, inforamses%ncpu
       Read(12,'(A13,I11)') dumchar, inforamses%ndim
       Read(12,'(A13,I11)') dumchar, inforamses%levelmin
       Read(12,'(A13,I11)') dumchar, inforamses%levelmax
       Read(12,'(A13,I11)') dumchar, inforamses%ngridmax
       Read(12,'(A13,I11)') dumchar, inforamses%nstep_coarse
       Read(12,*) 
       Read(12,'(A13,E24.15)') dumchar, inforamses%boxlen
       Read(12,'(A13,E24.15)') dumchar, inforamses%time
       Read(12,'(A13,E24.15)') dumchar, inforamses%aexp
       Read(12,'(A13,E24.15)') dumchar, inforamses%h0
       Read(12,'(A13,E24.15)') dumchar, inforamses%omega_m
       Read(12,'(A13,E24.15)') dumchar, inforamses%omega_l
       Read(12,'(A13,E24.15)') dumchar, inforamses%omega_k
       Read(12,'(A13,E24.15)') dumchar, inforamses%omega_b
       Read(12,'(A13,E24.15)') dumchar, inforamses%unit_l
       Read(12,'(A13,E24.15)') dumchar, inforamses%unit_d
       Read(12,'(A13,E24.15)') dumchar, inforamses%unit_t
       
       Close(12)
    End If

    Call Mpi_Bcast(inforamses, 1, Mpi_Type_Inforamses, 0, Mpi_Comm_World, ierr)


  End Subroutine readinforamses


  !=======================================================================
  Subroutine readinfoconepart(filename, procID, infocone)
    
    Implicit none

    Character(len=400), intent(in) :: filename
    Integer(kind=4), intent(in) :: procID
    Type(Type_infocone_part), intent(out) :: infocone

    ! Local variable
    Integer(kind=4) :: ierr
    Character(len=13) :: dumchar

    If(procID==0) Then
       
       Open(Unit=12, file=filename, status='old', iostat=ierr)
       If(ierr /= 0) Then
          Print *,'Error opening file ',trim(filename)
       End If
       Read(12,'(A13,I11)') dumchar, infocone%ncpu
       Read(12,'(A13,I11)') dumchar, infocone%nstride
       Read(12,'(A13,I11)') dumchar, infocone%nstep_coarse
       Read(12,'(A13,E24.15)') dumchar, infocone%aexp
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_x
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_y
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_z
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_rds
       Read(12,'(A13,I11)') dumchar, infocone%cone_id
       Read(12,'(A13,E24.15)') dumchar, infocone%cone_zlim
       Read(12,'(A13,E24.15)') dumchar, infocone%amax
       Read(12,'(A13,E24.15)') dumchar, infocone%amin
       Read(12,'(A13,E24.15)') dumchar, infocone%zmax
       Read(12,'(A13,E24.15)') dumchar, infocone%zmin
       Read(12,'(A13,E24.15)') dumchar, infocone%dmax
       Read(12,'(A13,E24.15)') dumchar, infocone%dmin
       Read(12,'(A13,E24.15)') dumchar, infocone%dtol
       Read(12,'(A13,I20)') dumchar, infocone%nglobalfile
       Read(12,'(A13,I20)') dumchar, infocone%npart
       Read(12,'(A13,I11)') dumchar, infocone%isfullsky
       Read(12,'(A13,E24.15)') dumchar, infocone%thetay
       Read(12,'(A13,E24.15)') dumchar, infocone%thetaz
       Read(12,'(A13,E24.15)') dumchar, infocone%theta
       Read(12,'(A13,E24.15)') dumchar, infocone%phi
       Read(12,'(A13,E24.15)') dumchar, infocone%aendconem2
       Read(12,'(A13,E24.15)') dumchar, infocone%aendconem1
       Read(12,'(A13,E24.15)') dumchar, infocone%aendcone
       Read(12,'(A13,E24.15)') dumchar, infocone%aexpold
       Read(12,'(A13,E24.15)') dumchar, infocone%aexp
       Read(12,'(A13,E24.15)') dumchar, infocone%zendconem2
       Read(12,'(A13,E24.15)') dumchar, infocone%zendconem1
       Read(12,'(A13,E24.15)') dumchar, infocone%zendcone
       Read(12,'(A13,E24.15)') dumchar, infocone%zexpold
       Read(12,'(A13,E24.15)') dumchar, infocone%zexp
       Read(12,'(A13,E24.15)') dumchar, infocone%dendconem2
       Read(12,'(A13,E24.15)') dumchar, infocone%dendconem1
       Read(12,'(A13,E24.15)') dumchar, infocone%dendcone
       Read(12,'(A13,E24.15)') dumchar, infocone%dexpold
       Read(12,'(A13,E24.15)') dumchar, infocone%dexp
       Read(12,'(A13,I11)') dumchar, infocone%future              
       Close(12)
    End If

    Call Mpi_Bcast(infocone, 1, Mpi_Type_infocone_part, 0, Mpi_Comm_World, ierr)

  End Subroutine readinfoconepart

  !=======================================================================
  Subroutine readinfoconegrav(filename, procID, infocone)
    
    Implicit none

    Character(len=400), intent(in) :: filename
    Integer(kind=4), intent(in) :: procID
    Type(Type_infocone_grav), intent(out) :: infocone

    ! Local variable
    Integer(kind=4) :: ierr
    Character(len=13) :: dumchar

    If(procID==0) Then
       
       Open(Unit=12, file=filename, status='old', iostat=ierr)
       If(ierr /= 0) Then
          Print *,'Error opening file ',trim(filename)
       End If
       Read(12,'(A13,I11)') dumchar, infocone%ncpu
       Read(12,'(A13,I11)') dumchar, infocone%nstride
       Read(12,'(A13,I11)') dumchar, infocone%nstep_coarse
       Read(12,'(A13,E24.15)') dumchar, infocone%aexp
       Read(12,'(A13,I11)') dumchar, infocone%nlevel
       Read(12,'(A13,I11)') dumchar, infocone%levelmin
       Read(12,'(A13,I11)') dumchar, infocone%levelmax
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_x
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_y
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_z
       Read(12,'(A13,E24.15)') dumchar, infocone%observer_rds
       Read(12,'(A13,I11)') dumchar, infocone%cone_id
       Read(12,'(A13,E24.15)') dumchar, infocone%cone_zlim
       Read(12,'(A13,E24.15)') dumchar, infocone%amax
       Read(12,'(A13,E24.15)') dumchar, infocone%amin
       Read(12,'(A13,E24.15)') dumchar, infocone%zmax
       Read(12,'(A13,E24.15)') dumchar, infocone%zmin
       Read(12,'(A13,E24.15)') dumchar, infocone%dmax
       Read(12,'(A13,E24.15)') dumchar, infocone%dmin
       Read(12,'(A13,E24.15)') dumchar, infocone%dtol
       Read(12,'(A13,I20)') dumchar, infocone%nglobalfile
       Read(12,'(A13,I20)') dumchar, infocone%nglobalcell
       Read(12,'(A13,I11)') dumchar, infocone%isfullsky
       Read(12,'(A13,E24.15)') dumchar, infocone%thetay
       Read(12,'(A13,E24.15)') dumchar, infocone%thetaz
       Read(12,'(A13,E24.15)') dumchar, infocone%theta
       Read(12,'(A13,E24.15)') dumchar, infocone%phi
       Read(12,'(A13,E24.15)') dumchar, infocone%aendconem2
       Read(12,'(A13,E24.15)') dumchar, infocone%aendconem1
       Read(12,'(A13,E24.15)') dumchar, infocone%aendcone
!!$       Read(12,'(A13,E24.15)') dumchar, infocone%aexpold
!!$       Read(12,'(A13,E24.15)') dumchar, infocone%aexp
       Read(12,'(A13,E24.15)') dumchar, infocone%zendconem2
       Read(12,'(A13,E24.15)') dumchar, infocone%zendconem1
       Read(12,'(A13,E24.15)') dumchar, infocone%zendcone
!!$       Read(12,'(A13,E24.15)') dumchar, infocone%zexpold
!!$       Read(12,'(A13,E24.15)') dumchar, infocone%zexp
       Read(12,'(A13,E24.15)') dumchar, infocone%dendconem2
       Read(12,'(A13,E24.15)') dumchar, infocone%dendconem1
       Read(12,'(A13,E24.15)') dumchar, infocone%dendcone
!!$       Read(12,'(A13,E24.15)') dumchar, infocone%dexpold
!!$       Read(12,'(A13,E24.15)') dumchar, infocone%dexp
       Read(12,'(A13,I11)') dumchar, infocone%future
       Close(12)
    End If

    Call Mpi_Bcast(infocone, 1, Mpi_Type_infocone_grav, 0, Mpi_Comm_World, ierr)

  End Subroutine readinfoconegrav

End Module modreadinfo
