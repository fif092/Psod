!==============================================================================
! Project: pFoF
! File: common/src/modhalo.f90
! Copyright Fabrice Roy and Vincent Bouillot (2011)
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
!! This file contains subroutines used to gather haloes on the differents processes and compute some observables.

!> This module contains subroutines used to gather haloes on the differents processes and compute some observables.
!>
!> Authors: F. Roy, V. Bouillot

Module modhalo

  Use modconstant

  Implicit none

  Integer(kind=4) :: haloNB        !< local number of haloes
  Integer(kind=4) :: haloNB_all    !< global number of haloes
  Integer(kind=4) :: halopartNB    !< number of particles belonging to one of the local haloes
  Integer(kind=4) :: final_local_npart  !< number of particles after the particles have been exchanged to gather them by haloes
  Integer(kind=4),   dimension(:), allocatable :: haloMass       !< mass of the haloes
  Integer(kind=PRI), dimension(:), allocatable :: haloID         !< ID of the haloes
  Integer(kind=PRI), dimension(:), allocatable :: halopartID     !< ID of the particles belonging to a halo
  Integer(kind=PRI), dimension(:), allocatable :: halopartramsesID     !< RAMSES ID of the particles belonging to a halo, used only for lightcone halo
  Integer(kind=4),   dimension(:), allocatable :: haloSubHaloNB  !< number of subhalo for each halo

  Real(kind=4), dimension(:,:), allocatable :: halopartPos    !< position of the particles belonging to one of the local haloes
  Real(kind=4), dimension(:,:), allocatable :: halopartVel    !< velocity of the particles belonging to one of the local haloes
  Real(kind=4), dimension(:), allocatable :: halopartPot    !< potential of the particles belonging to one of the local haloes
  Real(kind=4), dimension(:,:), allocatable :: halopartFor     !< force on the particles belonging to one of the local haloes
  Real(kind=8), dimension(:,:), allocatable :: halocomPos     !< position of the center of mass of each halo
  Real(kind=8), dimension(:,:), allocatable :: halocomVel     !< velocity of the center of mass of each halo
  Real(kind=8), dimension(:),   allocatable :: haloRadius     !< radius of each halo

  Private

  Public :: gatherhaloes, &
       selecthaloes, &
       computecom, &
       computeradius, &
       halopartPos, &
       halopartVel, &
       halopartPot, &
       halopartFor, &
       halocomPos, & 
       halocomVel, & 
       haloRadius, & 
       haloSubHaloNB, &! halopartramsesID 
       halopartID, & 
       haloID, & 
       haloMass, &
       haloNB, &
       haloNB_all, &
       halopartNB, &
       halopartramsesID, &
       final_local_npart

Contains

  !======================================================================

  Subroutine count_particles_number(local_npart, structure_id, values_count, particles_count)

    Use modmpicommons
    
    Implicit none
    
    Integer(kind=4), intent(in) :: local_npart
    Integer(kind=PRI), intent(in), dimension(:) :: structure_id
    Integer(kind=PRI), intent(inout), allocatable, dimension(:,:) :: particles_count

    Integer(kind=4) :: ip
    Integer(kind=PRI) :: current_id
    Integer(kind=4) :: current_value
    Integer(kind=4) :: values_count
    Integer(kind=4) :: unit
    
    current_id = structure_id(1)
    values_count = 1
    Do ip = 2, local_npart
       If(structure_id(ip) /= current_id) Then
          values_count = values_count + 1
          current_id = structure_id(ip)
       End If
    End Do

    Allocate(particles_count(2,values_count))

    current_value = 1
    particles_count(1,1) = structure_id(1)
    particles_count(2,1) = 1
    Do ip = 2, local_npart
       If(structure_id(ip) == particles_count(1,current_value)) Then
          particles_count(2,current_value) = particles_count(2,current_value) + 1
       Else
          current_value = current_value + 1
          particles_count(1,current_value) = structure_id(ip)
          particles_count(2,current_value) = 1          
       End If
    End Do

    
  End Subroutine count_particles_number
  
  !======================================================================
  
  Subroutine merge_particles_number(local_nvalues,particles_count)

    Use mpi
    Use modmpicommons
    
    Implicit none

    Integer(kind=4), intent(inout) :: local_nvalues
    Integer(kind=PRI), allocatable, dimension(:,:), intent(inout) :: particles_count
    Integer(kind=PRI), allocatable, dimension(:,:) :: merged_particles_count
    Integer(kind=PRI), allocatable, dimension(:,:) :: dest_particles_count

    Integer(kind=4) :: dest, right, left, dest_nvalues, merged_nvalues
    Integer(kind=4) :: loop, loop_number
    Logical(kind=4) :: even_pid
    Integer(kind=4) :: dest_ip, local_ip, merged_ip, first_ip, last_ip

    Integer(kind=4) :: mpierr
    Integer(kind=4), dimension(Mpi_Status_Size) :: mpistat

    Integer(kind=4) :: ih
    Integer(kind=4) :: local_halo_number
    Integer(kind=PRI), allocatable, dimension(:,:) :: local_halo_size_array
    
    Integer(kind=4) :: unit
    
    ! pid of the process right and left to current one
    right = mod(procID+1,procNB)
    left = mod(procID-1+procNB,procNB)

    ! first process has no left neighboor
    ! and last process has no right neighboor
    If(procID == 0) left = Mpi_Proc_Null
    If(procID == procNB-1) right = Mpi_Proc_Null

    ! is current process even?
    even_pid = (mod(procID,2)==0)
    
    ! there will be procNB/2 merging and 2 merging per loop => procNB/4 loops +1 for safety
    loop_number = procNB / 2 + 1

    
    Do loop = 1, loop_number

       ! merge 0 with 1, 2 with 3, 4 with 5, etc...
       If(even_pid) Then
          dest = right ! even process communicates with right neighboor
       Else
          dest = left  ! odd process communicates with left neighboor
       End If

       dest_nvalues = 0
       ! exchange of number of values in particles number per halo ID
       Call Mpi_Sendrecv(local_nvalues, 1, Mpi_Integer, dest, 0, &
            dest_nvalues, 1, Mpi_Integer, dest, 0, Mpi_Comm_World, mpistat, mpierr)

       ! allocate array for particles number per halo ID from neighbour
       If(dest_nvalues /= 0) Then
          Allocate(dest_particles_count(2,dest_nvalues))
          ! allocate array for merged particles number per halo ID
          merged_nvalues = (local_nvalues + dest_nvalues)
          Allocate(merged_particles_count(2,merged_nvalues))
          
          ! exchange particles number per halo ID 
          Call Mpi_Sendrecv(particles_count, 2*local_nvalues, Mpi_Pri, dest, 0, &
               dest_particles_count, 2*dest_nvalues, Mpi_Pri, dest, 0, Mpi_Comm_World, mpistat, mpierr)

          ! loop until we reach the end of both our and neighboor's array of particles number per halo ID
          local_ip = 1
          dest_ip = 1
          merged_ip = 1
          Do
             ! we hit the end of both array: end of the loop
             If(local_ip == local_nvalues .and. dest_ip == dest_nvalues) Exit
             ! we hit the end of local array: we copy the end of neighboor's array
             If(local_ip == local_nvalues) Then
                Print *,procID,' LOCAL TERMINE A MERGED=',merged_ip
                Do 
                   merged_particles_count(1,merged_ip) = dest_particles_count(1,dest_ip)
                   merged_particles_count(2,merged_ip) = dest_particles_count(2,dest_ip)
                   merged_ip = merged_ip + 1
                   dest_ip = dest_ip + 1
                   If(dest_ip > dest_nvalues) Exit
                End Do
                Exit
             End If
             ! we hit the end of neighboor's array: we copy the end of local array
             If(dest_ip == dest_nvalues) Then
                Print *,procID,' DEST TERMINE A MERGED=',merged_ip
                Do 
                   merged_particles_count(1,merged_ip) = particles_count(1,local_ip)
                   merged_particles_count(2,merged_ip) = particles_count(2,local_ip)
                   merged_ip = merged_ip + 1
                   local_ip = local_ip + 1
                   If(local_ip > local_nvalues) Exit
                End Do
                Exit
             End If
             ! we add next element from local array because the haloID is lower than neighboor's one
             If(particles_count(1,local_ip) < dest_particles_count(1,dest_ip)) Then
                merged_particles_count(1,merged_ip) = particles_count(1,local_ip) 
                merged_particles_count(2,merged_ip) = particles_count(2,local_ip)
                merged_ip = merged_ip + 1
                local_ip = local_ip + 1
                Cycle
                ! local and neighboor's halo ID are the same: same halo, we add particles numbers
             Else If(particles_count(1,local_ip) == dest_particles_count(1,dest_ip)) Then
                merged_particles_count(1,merged_ip) = particles_count(1,local_ip) 
                merged_particles_count(2,merged_ip) = particles_count(2,local_ip) + &
                     dest_particles_count(2,dest_ip)
                merged_ip = merged_ip + 1
                local_ip = local_ip + 1
                dest_ip = dest_ip + 1
                Cycle
                ! we add next element from neighboor's array because the haloID is lower than local one
             Else If(particles_count(1,local_ip) > dest_particles_count(1,dest_ip)) Then
                merged_particles_count(1,merged_ip) = dest_particles_count(1,dest_ip) 
                merged_particles_count(2,merged_ip) = dest_particles_count(2,dest_ip)
                merged_ip = merged_ip + 1
                dest_ip = dest_ip + 1
                Cycle
             End If
          End Do
          
          ! decread merge_ip by one (correction from the loop)
          merged_ip = merged_ip - 1
       
          ! even: keep left part of the merged count array
          If(even_pid) Then
             local_nvalues = merged_ip/2
             If(mod(merged_ip,2)==1) local_nvalues = local_nvalues+1
             first_ip = 1
             last_ip = local_nvalues
          Else ! odd: keep right part of the merged count array
             local_nvalues = merged_ip/2
             first_ip = local_nvalues+1
             If(mod(merged_ip,2)==1) first_ip = first_ip+1
             last_ip = first_ip + local_nvalues - 1
          End If

          Deallocate(particles_count)
          Allocate(particles_count(2,local_nvalues))
          particles_count(:,:) = merged_particles_count(:,first_ip:last_ip)
          Deallocate(merged_particles_count)
       
          Deallocate(dest_particles_count)

       End If


       ! 1 with 2, 3 with 4, etc...
       If(even_pid) Then
          dest = left
       Else
          dest = right
       End If

       dest_nvalues = 0
       ! exchange of number of values in particles number per halo ID
       Call Mpi_Sendrecv(local_nvalues, 1, Mpi_Integer, dest, 0, &
            dest_nvalues, 1, Mpi_Integer, dest, 0, Mpi_Comm_World, mpistat, mpierr)

       If(dest_nvalues /= 0) Then
          ! allocate array for particles number per halo ID from neighbour
          Allocate(dest_particles_count(2,dest_nvalues))
          ! allocate array for merged particles number per halo ID
          merged_nvalues = (local_nvalues + dest_nvalues)
          Allocate(merged_particles_count(2,merged_nvalues))
          
          ! exchange particles number per halo ID 
          Call Mpi_Sendrecv(particles_count, 2*local_nvalues, Mpi_Pri, dest, 0, &
               dest_particles_count, 2*dest_nvalues, Mpi_Pri, dest, 0, Mpi_Comm_World, mpistat, mpierr)
          
          ! loop until we reach the end of both our and neighboor's array of particles number per halo ID
          local_ip = 1
          dest_ip = 1
          merged_ip = 1
          Do
             ! we hit the end of both array: end of the loop
             If(local_ip == local_nvalues .and. dest_ip == dest_nvalues) Exit
             ! we hit the end of local array: we copy the end of neighboor's array
             If(local_ip == local_nvalues) Then
                Do 
                   merged_particles_count(1,merged_ip) = dest_particles_count(1,dest_ip)
                   merged_particles_count(2,merged_ip) = dest_particles_count(2,dest_ip)
                   merged_ip = merged_ip + 1
                   dest_ip = dest_ip + 1
                   If(dest_ip > dest_nvalues) Exit
                End Do
                Exit
             End If
             ! we hit the end of neighboor's array: we copy the end of local array
             If(dest_ip == dest_nvalues) Then
                Do 
                   merged_particles_count(1,merged_ip) = particles_count(1,local_ip)
                   merged_particles_count(2,merged_ip) = particles_count(2,local_ip)
                   merged_ip = merged_ip + 1
                   local_ip = local_ip + 1
                   If(local_ip > local_nvalues) Exit
                End Do
                Exit
             End If
             ! we add next element from local array because the haloID is lower than neighboor's one
             If(particles_count(1,local_ip) < dest_particles_count(1,dest_ip)) Then
                merged_particles_count(1,merged_ip) = particles_count(1,local_ip) 
                merged_particles_count(2,merged_ip) = particles_count(2,local_ip)
                merged_ip = merged_ip + 1
                local_ip = local_ip + 1
                Cycle
                ! local and neighboor's halo ID are the same: same halo, we add particles numbers
             Else If(particles_count(1,local_ip) == dest_particles_count(1,dest_ip)) Then
                merged_particles_count(1,merged_ip) = particles_count(1,local_ip) 
                merged_particles_count(2,merged_ip) = particles_count(2,local_ip) + &
                     dest_particles_count(2,dest_ip)
                merged_ip = merged_ip + 1
                local_ip = local_ip + 1
                dest_ip = dest_ip + 1
                Cycle
                ! we add next element from neighboor's array because the haloID is lower than local one
             Else If(particles_count(1,local_ip) > dest_particles_count(1,dest_ip)) Then
                merged_particles_count(1,merged_ip) = dest_particles_count(1,dest_ip) 
                merged_particles_count(2,merged_ip) = dest_particles_count(2,dest_ip)
                merged_ip = merged_ip + 1
                dest_ip = dest_ip + 1
                Cycle
             End If
          End Do
          
          ! decread merge_ip by one (correction from the loop)
          merged_ip = merged_ip - 1
          
          ! even: keep right part of the merged count array
          If(even_pid) Then
             local_nvalues = merged_ip/2
             first_ip = local_nvalues+1
             If(mod(merged_ip,2)==1) first_ip = first_ip+1
             last_ip = first_ip + local_nvalues - 1
          Else ! odd: keep left part of the merged count array
             local_nvalues = merged_ip/2
             If(mod(merged_ip,2)==1) local_nvalues = local_nvalues+1
             first_ip = 1
             last_ip = local_nvalues
          End If
          
          Deallocate(particles_count)
          Allocate(particles_count(2,local_nvalues))
          particles_count(:,:) = merged_particles_count(:,first_ip:last_ip)
          Deallocate(merged_particles_count)
          
          Deallocate(dest_particles_count)
       End If
       
    End Do
    
    ! count the number of halo with number of particles > Mmin
!    local_halo_number = 0
!    Do local_ip = 1, local_nvalues
!       If(particles_count(2,local_ip) >= param%mmin) Then
!          local_halo_number = local_halo_number + 1
!       End If
!    End Do

    ! allocate array to store the number of particles per halo and the corresponding haloID
!    Allocate(local_halo_size_array(2,local_halo_number))

    ! we only keep halo ID and particles number for halo with particles number > Mmin
!    ih = 1
!    Do local_ip = 1, local_nvalues
!       If(particles_count(2,local_ip) >= param%mmin) Then
!          local_halo_size_array(1,ih) = particles_count(1,local_ip)
!          local_halo_size_array(2,ih) = particles_count(2,local_ip)
!          ih = ih + 1
!       End If
!    End Do

   
  End Subroutine merge_particles_number
  
  !======================================================================
  Subroutine gatherhaloLB()
    
    Use modvarcommons
    Use modsort

    Implicit none

    Integer(kind=4) :: local_nvalues
    Integer(kind=PRI), allocatable, dimension(:,:) :: particles_count
    
    ! sort the particles according structure_id
    Call heapsort(local_npart, structure_id, position, velocity, field, potential, pfof_id) !, ramses_id)

    Call count_particles_number(local_npart, structure_id, local_nvalues,particles_count)

    Call merge_particles_number(local_nvalues,particles_count)
    
  End Subroutine gatherhaloLB
  
  
  !=======================================================================
  !> Exchange the particles so that particles belonging to one halo are gathered on the same process.
  Subroutine gatherhaloes(mpicomm, param)

    Use mpi
    Use modvarcommons
    Use modmpicommons
    Use modtiming

    Implicit None

    Integer, intent(in) :: mpicomm   !< MPI communicator used for the communications
    Class(Type_parameter_pfof), intent(in) :: param

    Integer(kind=PRI), dimension(:), allocatable :: strSend ! structure ID array, and tmp array used for comm
    Integer(kind=PRI), dimension(:), allocatable :: idSend  ! particle ID array and tmp array used for comm
    Integer(kind=PRI), dimension(:), allocatable :: ramsesidSend
    Real(kind=4), dimension(:,:), allocatable :: posSend, velSend, forSend
    Real(kind=4), dimension(:), allocatable :: potSend
    Integer(kind=4) :: strPID
    Integer(kind=4), dimension(:),allocatable :: strPIDvec, strPIDvecloc  ! array of particle nb by process after
    Integer(kind=4) :: sendID, recvID   ! ID for process to send to and to receive from
    Integer(kind=4) :: i, ind, iproc
    Integer(kind=4) :: allocStat
    Integer(kind=4) :: mpistat(MPI_STATUS_SIZE)   ! MPI comm. status
    Integer(kind=4) :: nbrec, nbsend    ! nb of elements to recv and to send
    Integer(kind=PRI) :: strNBbeg, strNBend
    Integer(kind=PRI) :: smin,smax
    Integer(kind=PRI) :: tmpdi
    Integer(kind=4) :: NPparProc
    Integer(kind=4) :: recvpoint
    Integer(kind=4) :: mpierr
    Integer(kind=4) :: mpireqs1, mpireqs2, mpireqs3, mpireqs4, mpireqs5, mpireqs6, mpireqs7
    Integer(kind=4) :: mpireqr1, mpireqr2, mpireqr3, mpireqr4, mpireqr5, mpireqr6, mpireqr7

    Logical(kind=4) :: do_read_ramses_part_id

    !!Call gatherhaloLB()
    
    Select Type (param)
    Type is (Type_parameter_pfof_snap)
       do_read_ramses_part_id = .false.
    Type is (Type_parameter_pfof_cone)
       do_read_ramses_part_id = param%do_read_ramses_part_id
    End Select

    tmpdi = global_npart/procNB
    NPparProc = int(tmpdi, kind=4)
    If(mod(global_npart,int(procNB,kind=8)) /= 0) Then
       NPparProc = NPparProc + 1
    End If
       
    strNBbeg = int(NPparProc,kind=8) * int(procID,kind=8) + 1    ! id min de halo sur le process courant
    strNBend = int(NPparProc,kind=8) * int((procID + 1),kind=8)  ! id max de halo sur le process courant
    
    Allocate(strPIDvecloc(procNB), strPIDvec(procNB))
    strPIDvecloc = 0  ! nb de particules par process
    strPIDvec = 0

    Do i = 1,local_npart
       strPID = int((structure_id(i)-1) / NPparProc, kind=4) + 1  ! id du process ou se trouvera la particule
       strPIDvecloc(strPID) = strPIDvecloc(strPID) + 1         ! on incremente le nb de particules sur ce process
    End Do
    
    Call Mpi_Allreduce(strPIDvecloc,strPIDvec,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

    final_local_npart = strPIDvec(procID+1)             ! nb de particules sur le process courant apres redistribution selon les halos
    Allocate(fposition(3,final_local_npart),STAT=allocStat)      
    If(allocStat > 0) Call EmergencyStop('Allocate failed for fposition in gatherhaloes',2)
    Allocate(fvelocity(3,final_local_npart),STAT=allocStat)
    If(allocStat > 0) Call EmergencyStop('Allocate failed for fvelocity in gatherhaloes',2)
    Allocate(fstructure_id(final_local_npart),STAT=allocStat)
    If(allocStat > 0) Call EmergencyStop('Allocate failed for fstructure_id in gatherhaloes',2)
    Allocate(fpfof_id(final_local_npart),STAT=allocStat)
    If(allocStat > 0) Call EmergencyStop('Allocate failed for fpfof_id in gatherhaloes',2)
    If(param%do_read_potential) Then
       Allocate(fpotential(final_local_npart),STAT=allocStat)
       If(allocStat > 0) Call EmergencyStop('Allocate failed for fpotential in gatherhaloes',2)
    End If
    If(param%do_read_gravitational_field) Then
       Allocate(ffield(3,final_local_npart),STAT=allocStat)
       If(allocStat > 0) Call EmergencyStop('Allocate failed for ffield in gatherhaloes',2)
    End If
    If(do_read_ramses_part_id) Then
       Allocate(framses_id(final_local_npart),STAT=allocStat)
       If(allocStat > 0) Call EmergencyStop('Allocate failed for framses_id in gatherhaloes',2)
    End If

    recvpoint = 1
    procMasse : Do iproc = 1,procNB-1
       If(procID==0) Write(*,'(A,I5,A)') 'Gatherhalo permutation ',iproc,'...'

       sendID = mod(procID + iproc,procNB)
       recvID = mod(procID + procNB - iproc,procNB)
       nbsend = strPIDvecloc(sendID+1)

       Call Mpi_ISend(nbsend,1,Mpi_Integer,sendID,sendID,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(nbrec, 1,Mpi_Integer,recvID,procID,mpicomm,mpireqr1,mpierr)

       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       If(nbsend /= 0) Then
          Allocate(strSend(nbsend),posSend(3,nbsend),velSend(3,nbsend),idSend(nbsend))
          If(param%do_read_potential) Allocate(potSend(nbsend))
          If(param%do_read_gravitational_field) Allocate(forSend(3,nbsend))
          If(do_read_ramses_part_id) Allocate(ramsesidSend(nbsend))

          smin = int(NPparProc,kind=8) * int(sendID,kind=8) + 1
          smax = int(NPparProc,kind=8) * int((sendID+1),kind=8)

          ind=1
          Do i=1, local_npart
             If(structure_id(i)>= smin .and. structure_id(i)<= smax) Then
                strSend(ind) = structure_id(i)
                posSend(:,ind) = position(:,i)
                velSend(:,ind) = velocity(:,i)
                If(param%do_read_potential) potSend(ind) = potential(i)
                If(param%do_read_gravitational_field) forSend(:,ind) = field(:,i)
                If(do_read_ramses_part_id) ramsesidSend(ind) = ramses_id(i)
                idSend(ind) = pfof_id(i)
                ind = ind+1
             End If
          End Do

          If(ind /= nbsend +1 ) Then
             write(*,*)ind,nbsend,sendID,smin,smax,strPIDvec
             Call EmergencyStop('Error  1 while sharing structures for output.',2)
          End If

          Call Mpi_Isend(strSend,nbsend,  MPI_PRI, sendID,1,mpicomm,mpireqs1,mpierr)
          Call Mpi_Isend(posSend,3*nbsend,Mpi_Real,sendID,2,mpicomm,mpireqs2,mpierr)
          Call Mpi_Isend(velSend,3*nbsend,Mpi_Real,sendID,3,mpicomm,mpireqs3,mpierr)
          Call Mpi_Isend(idSend, nbsend,  MPI_PRI, sendID,4,mpicomm,mpireqs4,mpierr)
          If(param%do_read_potential) Call Mpi_Isend(potSend, nbsend, Mpi_Real, sendID, 5,mpicomm,mpireqs5,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_Isend(forSend,3*nbsend,Mpi_Real,sendID,6,mpicomm,mpireqs6,mpierr)
          If(do_read_ramses_part_id) Call Mpi_Isend(ramsesidSend,nbsend,MPI_PRI,&
               sendID,7,mpicomm,mpireqs7,mpierr)

       End If

       If(nbrec /= 0) Then
          Call Mpi_IRecv(fstructure_id(recvpoint), nbrec,  MPI_PRI, recvID,1,mpicomm,mpireqr1,mpierr)
          Call Mpi_IRecv(fposition(1,recvpoint),3*nbrec,Mpi_Real,recvID,2,mpicomm,mpireqr2,mpierr)
          Call Mpi_IRecv(fvelocity(1,recvpoint),3*nbrec,Mpi_Real,recvID,3,mpicomm,mpireqr3,mpierr)
          Call Mpi_IRecv(fpfof_id(recvpoint), nbrec,  MPI_PRI, recvID,4,mpicomm,mpireqr4,mpierr)
          If(param%do_read_potential) Call Mpi_IRecv(fpotential(recvpoint), nbrec,  &
               Mpi_Real, recvID,5,mpicomm,mpireqr5,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_IRecv(ffield(1,recvpoint),3*nbrec, &
               Mpi_Real,recvID,6,mpicomm,mpireqr6,mpierr)
          If(do_read_ramses_part_id) Call Mpi_Irecv(framses_id(recvpoint), nbrec, MPI_PRI, &
               recvID, 7, mpicomm, mpireqr7, mpierr)
          recvpoint=recvpoint+nbrec
       End If


       If(nbsend/=0) Then
          Call Mpi_Wait(mpireqs1,mpistat,mpierr)
          Call Mpi_Wait(mpireqs2,mpistat,mpierr)
          Call Mpi_Wait(mpireqs3,mpistat,mpierr)
          Call Mpi_Wait(mpireqs4,mpistat,mpierr)
          Deallocate(strSend,posSend,idSend,velSend)

          If(param%do_read_potential) Then
             Call Mpi_Wait(mpireqs5,mpistat,mpierr)
             Deallocate(potSend)
          End If
          If(param%do_read_gravitational_field) Then
             Call Mpi_Wait(mpireqs6,mpistat,mpierr)
             Deallocate(forSend)
          End If
          If(do_read_ramses_part_id) Then
             Call Mpi_Wait(mpireqs7, mpistat, mpierr)
             Deallocate(ramsesidSend)
          End If
          
       End If
       If(nbrec/=0) Then
          Call Mpi_Wait(mpireqr1,mpistat,mpierr)
          Call Mpi_Wait(mpireqr2,mpistat,mpierr)
          Call Mpi_Wait(mpireqr3,mpistat,mpierr)
          Call Mpi_Wait(mpireqr4,mpistat,mpierr)
          If(param%do_read_potential) Call Mpi_Wait(mpireqr5,mpistat,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_Wait(mpireqr6,mpistat,mpierr)
          If(do_read_ramses_part_id) Call Mpi_Wait(mpireqr7, mpistat, mpierr)
       End If

    End Do procMasse


    ind= recvpoint
    Do i=1, local_npart
       If(structure_id(i)>= strNBbeg .and. structure_id(i)<= strNBend) Then
          fstructure_id(recvpoint)  = structure_id(i)
          fposition(:,recvpoint) = position(:,i)
          fvelocity(:,recvpoint) = velocity(:,i)
          fpfof_id(recvpoint)  = pfof_id(i)
          If(param%do_read_potential) fpotential(recvpoint) = potential(i)
          If(param%do_read_gravitational_field) ffield(:,recvpoint) = field(:,i)
          If(do_read_ramses_part_id) framses_id(recvpoint) = ramses_id(i)
          recvpoint = recvpoint+1
       End If
    End Do

    Deallocate(structure_id)
    Deallocate(strPIDvecloc,strPIDvec)
    Deallocate(position, velocity, pfof_id)
    If(param%do_read_potential) Deallocate(potential)
    If(param%do_read_gravitational_field) Deallocate(field)
    If(do_read_ramses_part_id) Deallocate(ramses_id)

    If( recvpoint/= final_local_npart +1 ) Then
       Print *, procID,' recvpoint=',recvpoint
       Print *, procID,' final_local_npart=',final_local_npart
       Call EmergencyStop('Error while gathering haloes in gatherhaloes.',2)
    End If


  End Subroutine gatherhaloes


  ! ======================================================================
  ! Select haloes whose mass is >= Mmin
  Subroutine selecthaloes(param)

#ifdef DEBUG
    Use mpi
#endif
    Use modconstant
    Use modvarcommons
    Use modmpicommons, only : procID, EmergencyStop

    Implicit none

    Class(Type_parameter_pfof), intent(in) :: param

    Integer(kind=8) :: hidmin, hidmax
    Integer(kind=4) :: nbhid
    Integer(kind=4) :: hindex
    Integer(kind=4) :: li, lp, fi, fp, h, i
    Integer(kind=4), dimension(:), allocatable :: halomasstmp

    Logical(kind=4) :: do_read_ramses_part_id


    Select Type (param)
    Type is (Type_parameter_pfof_snap)
       do_read_ramses_part_id = .false.
    Type is (Type_parameter_pfof_cone)
       do_read_ramses_part_id = param%do_read_ramses_part_id
    End Select

#ifdef DEBUG
    Print *, 'Process ',procID, ' enters selecthaloes'
#endif
    
    If(final_local_npart /= 0) Then
       hidmin = fstructure_id(1)
       hidmax = fstructure_id(final_local_npart)
    Else
       hidmin=0
       hidmax=-1
    End If
    nbhid = int(hidmax - hidmin + 1, kind=4)

#ifdef DEBUG
    Print *, 'hidmax=',hidmax,' ; hidmin=',hidmin,' ; nbhid=',nbhid, ' ; final_local_npart=',final_local_npart 
#endif

    Allocate(halomasstmp(nbhid))

    halomasstmp = 0

    Do i=1, final_local_npart
       hindex = int(fstructure_id(i) - hidmin + 1,kind=4)
       halomasstmp(hindex) = halomasstmp(hindex) + 1
    End Do

    haloNB = 0
    halopartNB = 0

    ! Compute total nb of particles in halo with M >= Mmin and nb of halos with M >= Mmin
    Do i = 1, nbhid
       If(halomasstmp(i) >= param%mmin) Then
          haloNB = haloNB + 1
          halopartNB = halopartNB + halomasstmp(i)
       End If
    End Do


    ! Keep positions, velocities and id for particles in halo with M >= Mmin, and potential if requested
    Allocate(halopartPos(3,halopartNB))
    Allocate(halopartVel(3,halopartNB))
    Allocate(halopartID(halopartNB))
    If(param%do_read_potential) Allocate(halopartPot(halopartNB))
    If(param%do_read_gravitational_field) Allocate(halopartFor(3,halopartNB))
    If(do_read_ramses_part_id) Allocate(halopartramsesid(halopartNB))
    ! Keep mass and id for halos with M >= Mmin
    If(haloNB==0) Then
       Allocate(haloMass(1))
       Allocate(haloID(1))
    Else
       Allocate(haloMass(haloNB))
       Allocate(haloID(haloNB))
    End If
    ! Sub-halo detection is not implemented yet
    ! We allocate halosubhaloNB with a size=1
    Allocate(halosubhaloNB(1))


    fp = 1
    lp = 0
    fi = 1
    li = 0
    h = 1
    Do i = 1, nbhid
       If(halomasstmp(i) >= param%mmin) Then
          lp = fp + halomasstmp(i) - 1
          li = fi + halomasstmp(i) - 1
          halopartPos(:,fp:lp) = fposition(:,fi:li)
          halopartVel(:,fp:lp) = fvelocity(:,fi:li)
          If(param%do_read_potential) halopartPot(fp:lp) = fpotential(fi:li)
          If(param%do_read_gravitational_field) halopartFor(:,fp:lp) = ffield(:,fi:li)
          If(do_read_ramses_part_id) halopartramsesid(fp:lp) = framses_id(fi:li)
          halopartID(fp:lp) = fpfof_id(fi:li)
          haloMass(h) = halomasstmp(i)
          haloID(h) = fstructure_id(fi)
          fp = lp + 1
          fi = li + 1
          h = h + 1
       Else
          li = fi + halomasstmp(i) - 1
          fi = li + 1
       End If
    End Do
    If(li /= final_local_npart) Then
       Print *, 'Error when keeping particles positions in halo with M >= Mmin on process ', procID
    End If
    If(lp /= halopartNB) Then
       Print *, 'Error when keeping particles positions in halo with M >= Mmin on process ', procID
    End If


    If(Allocated(fposition)) Deallocate(fposition)
    If(Allocated(fvelocity)) Deallocate(fvelocity)
    If(Allocated(fpfof_id)) Deallocate(fpfof_id)
    If(Allocated(halomasstmp)) Deallocate(halomasstmp)
    If(Allocated(fpotential)) Deallocate(fpotential)
    If(Allocated(ffield)) Deallocate(ffield)
    If(Allocated(framses_id)) Deallocate(framses_id)
  End Subroutine selecthaloes


  ! ======================================================================
  ! Computes the position and the velocity of the center of mass for each halo
  Subroutine computecom(periodic)

    Use mpi
    Implicit none

    Logical, intent(in) :: periodic

    Integer(kind=4) :: fi, li, h, i, j
    Integer(kind=4) :: halom
    Real(kind=8), dimension(3) :: delta
    Integer(kind=4) :: mpierr

    If(haloNB==0) Then
       Allocate(halocomPos(3,1))
       Allocate(halocomVel(3,1))
    Else
       Allocate(halocomPos(3,haloNB))
       Allocate(halocomVel(3,haloNB))
    End If

    halocomPos = 0.d0
    halocomVel = 0.d0

    fi = 1
    Do h = 1, haloNB
       halom = 1
       li = fi + haloMass(h) - 1
       halocomPos(:,h) = halopartPos(:,fi)
       halocomVel(:,h) = halopartVel(:,fi)
       Do i = fi+1, li
          delta(:) = halocomPos(:,h) / halom - halopartPos(:,i)
          If(periodic) Then
             Do j = 1, 3
                If(abs(delta(j)) > 0.5d0) Then
                   If(delta(j) > 0d0) halocomPos(j,h) = halocomPos(j,h) + 1.0d0
                   If(delta(j) < 0d0) halocomPos(j,h) = halocomPos(j,h) - 1.0d0
                End If
             End Do
          End If
          halom = halom + 1
          halocomPos(:,h) = halocomPos(:,h) + halopartPos(:,i)
          halocomVel(:,h) = halocomVel(:,h) + halopartVel(:,i)
       End Do
       fi = li + 1
       
    End Do

    Do h = 1, haloNB
       halocomPos(:,h) = halocomPos(:,h) / haloMass(h)
       halocomVel(:,h) = halocomVel(:,h) / haloMass(h)
       If(periodic) Then
          Do j = 1, 3
             If(halocomPos(j,h) > 1.0d0) halocomPos(j,h) = halocomPos(j,h) - 1.0d0
             If(halocomPos(j,h) < 0.0d0) halocomPos(j,h) = halocomPos(j,h) + 1.0d0
          End Do
       End If
    End Do
    
    Call Mpi_AllReduce(haloNB,haloNB_all,1,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

  End Subroutine computecom

  
  ! ======================================================================
  Subroutine computeradius(periodic)

    Implicit None
    
    Logical, intent(in) :: periodic

    Integer(kind=4) :: ib
    Integer(kind=4) :: id
    Integer(kind=4) :: ih
    Integer(kind=4) :: ip

    Real(kind=8) :: rmax2
    Real(kind=8) :: d2
    Real(kind=8), dimension(3) :: delta
    
    If(haloNB==0) Then
       Allocate(haloRadius(1))
    Else
       Allocate(haloRadius(haloNB))
    End If

    ib = 0
    Do ih = 1, haloNB
       rmax2 = 0.d0
       Do ip = ib+1, ib+haloMass(ih)
          delta(:) = halopartPos(:,ip) - haloComPos(:,ih)
          If(periodic) Then
             Do id = 1, 3
                If(abs(delta(id)) > 0.5d0) Then
                   delta(id) = 1.d0 - abs(delta(id))
                End If
             End Do
          End If
          d2 = delta(1)*delta(1) + delta(2)*delta(2) + delta(3)*delta(3) 
          If(d2 > rmax2) rmax2 = d2
       End Do
       ib = ib + haloMass(ih)
       haloRadius(ih) = sqrt(rmax2)
    End Do

    Deallocate(halopartPos, halopartVel)


  End Subroutine computeradius

  
End Module modhalo
