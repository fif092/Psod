!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modextract
  
  Use modvariables

  Implicit none
  
  Real(kind=4), dimension(3) :: center 
  Character(len=8) :: type                   ! sphere, cube, cuboid
  Real(kind=4), dimension(3) :: dimensions   ! (radius, 0,0) for a sphere, (a,0,0) for a cube, (dx,dy,dz) for a cuboid
  Character(len=400) :: filename
  Logical(kind=4) :: potential = .false.
  Logical(kind=4) :: force = .false.


Contains

  Subroutine extract_cuboid()
    
    Use modhdf5
    Implicit none

    Character(len=H5STRLEN) :: aname
    Character(len=H5STRLEN) :: dname
    Character(len=H5STRLEN) :: grname
    Character(len=8) :: charigroup
    Integer(kind=4) :: ifile
    Real(kind=4), dimension(:,:), allocatable :: boundaries
    Character(len=5) :: charifile
    Integer(kind=4) :: nfile
    Integer(kind=4) :: ngroup
    Logical(kind=4) :: fileexist
    Character(len=400) :: filehdf5
    Integer(kind=hid_t) :: file_id, gr_id, tmpid, dataid,metaid
    Integer(kind=4) :: groupsizem1

    Integer(kind=4) :: i,j,k
    Integer(kind=4) :: iloc, jloc, kloc
    Integer(kind=4) :: ind, indpart
    Integer(kind=4) :: ncglob
    Integer(kind=4) :: nc
    Integer(kind=4), dimension(:,:), allocatable :: globindexcube
    Integer(kind=4), dimension(6) :: globindexextract
    Integer(kind=4), dimension(3) :: dimcube2extract
    Integer(kind=4), dimension(:,:,:), allocatable :: coord2file
    Integer(kind=4) :: imin, imax, jmin, jmax, kmin, kmax
    Integer(kind=4), dimension(:,:), allocatable :: locindex2extract
    Integer(kind=4), dimension(:), allocatable :: ngroupperfile

    Integer(kind=4), dimension(:,:), allocatable :: npartpergroup
    Integer(kind=4) :: groupid, igroup
    Integer(kind=4) :: npart2extract
    Integer(kind=4) :: ngroup2extract

    Integer(kind=4) :: npart
    Integer(kind=4) :: potentialattr
    Integer(kind=4) :: forceattr

    ! If you use hdf5 earlier in your code you can comment this
    Call hdf5_init()

    filehdf5 = trim(filename)//'_00000.h5'

    ! Open the first hdf5 cube file and read some info
    Call hdf5_open_file(filehdf5, file_id)

    ! open the root group
    grname = 'metadata'
    Call hdf5_open_group(file_id,grname, metaid)

    aname='nfile'
    Call hdf5_read_attr(metaid,aname,nfile)

    aname='ngroup'
    Call hdf5_read_attr(metaid,aname,ngroup)

    aname='1/groupsize'
    Call hdf5_read_attr(metaid,aname,groupsizem1)

    Allocate(boundaries(6,nfile))
    aname = 'boundaries'
    Call hdf5_read_attr(metaid,aname,6,boundaries(:,1))

    grname='pfof_snap_parameters'
    Call hdf5_open_group(metaid, grname, tmpid)
    grname='input_parameters'
    Call hdf5_open_group(tmpid, grname, gr_id)

    aname='do_read_potential'
    Call hdf5_read_attr(gr_id,aname,potentialattr)
    If(potential .And. (potentialattr==0)) Then
       Print *,'Potential requested but not found in hdf5 file. Exit.'
       Stop 10
    End If

    aname='do_read_gravitational_field'
    Call hdf5_read_attr(gr_id,aname,forceattr)
    If(force .And. (forceattr==0)) Then
       Print *,'Force requested but not found in hdf5 file. Exit.'
       Stop 10
    End If

    Call hdf5_close_group(gr_id)
    Call hdf5_close_group(tmpid)
 
    Allocate(npartpergroup(ngroup,nfile))
    
    dname = 'npart_grp_array'
    Call hdf5_read_data(metaid,dname,ngroup,npartpergroup(:,1))
    
    Call hdf5_close_group(metaid)
    Call hdf5_close_file(file_id)

    Do ifile=1,nfile-1
       ! we open the files and read the boundaries (we have to, in case the set of files is not complete)
       Write(charifile(1:5),'(I5.5)') ifile
       filehdf5=trim(filename)//'_'//charifile//'.h5'
       Inquire(file=filehdf5, exist=fileexist)
       If(fileexist) Then
          Call hdf5_open_file(filehdf5,file_id)
          grname='metadata'
          Call hdf5_open_group(file_id,grname,metaid)
          aname = 'boundaries'
          Call hdf5_read_attr(metaid,aname,6,boundaries(:,ifile+1))

          dname = 'npart_grp_array'
          Call hdf5_read_data(metaid,dname,ngroup,npartpergroup(:,ifile+1))
          Call hdf5_close_group(metaid)
          Call hdf5_close_file(file_id)
       End If
    End Do

    Allocate(globindexcube(6,nfile))

    ! nb of "group cube" in each direction in the global domain : ngroup*nfile = ncglob**3
    ncglob = groupsizem1
    nc = int(ngroup**(1./3.),kind=4)

    globindexcube(:,:) = int( boundaries(:,:)*ncglob, kind=4 ) + 1
    ! the x,y,z-max are excluded, the index should be decreased by 1
    globindexcube(2,:) = globindexcube(2,:) - 1
    globindexcube(4,:) = globindexcube(4,:) - 1
    globindexcube(6,:) = globindexcube(6,:) - 1
    
    Do i=1,3
       globindexextract(2*i-1) = int((center(i) - 0.5*dimensions(i)) * ncglob, kind=4) + 1
       globindexextract(2*i)   = int((center(i) + 0.5*dimensions(i)) * ncglob, kind=4) + 1
       dimcube2extract(i) = globindexextract(2*i)-globindexextract(2*i-1)+1
    End Do
    
    Allocate(coord2file(ncglob, ncglob, ncglob))

    Do ifile=1,nfile
       imin = globindexcube(1,ifile)
       imax = globindexcube(2,ifile)
       jmin = globindexcube(3,ifile)
       jmax = globindexcube(4,ifile)
       kmin = globindexcube(5,ifile)
       kmax = globindexcube(6,ifile)
       ! warning: 1 <= ifile <= nfile
       coord2file(imin:imax,jmin:jmax,kmin:kmax) = ifile
    End Do

    ngroup2extract = dimcube2extract(1)*dimcube2extract(2)*dimcube2extract(3)

    ind=1
    Allocate(locindex2extract(2,ngroup2extract))
    Allocate(ngroupperfile(nfile))

    npart2extract = 0
    ngroupperfile = 0
    Do i=globindexextract(1), globindexextract(2)
       Do j=globindexextract(3), globindexextract(4)
          Do k=globindexextract(5), globindexextract(6)
             ifile = coord2file(i,j,k)
             ngroupperfile(ifile) = ngroupperfile(ifile)+1
             locindex2extract(2,ind) = ifile
             iloc = i - globindexcube(1,ifile) + 1
             jloc = j - globindexcube(3,ifile) + 1
             kloc = k - globindexcube(5,ifile) + 1
             locindex2extract(1,ind) = iloc + (jloc-1)*nc +(kloc-1)*nc*nc
             npart2extract = npart2extract + npartpergroup(locindex2extract(1,ind),ifile)
             ind = ind+1
          End Do
       End Do
    End Do
    
    Call quicksort(1,ngroup2extract,locindex2extract)

    Allocate(pos(3,npart2extract), vel(3,npart2extract), id(npart2extract))
    If(potential) Allocate(pot(npart2extract))

    ind = 0
    indpart = 1
    Do ifile=1,nfile
       If(ngroupperfile(ifile) /= 0) Then ! there is at least 1 group to read in this file
          Write(charifile(1:5),'(I5.5)') ifile-1
          filehdf5=trim(filename)//'_'//charifile//'.h5'
          Inquire(file=filehdf5, exist=fileexist)
          If(fileexist) Then
             Call hdf5_open_file(filehdf5,file_id)
             grname='data'
             Call hdf5_open_group(file_id,grname,dataid)
             ! loop over the ngroup to read in this file
             Do igroup=1,ngroupperfile(ifile)
                ind=ind+1
                groupid= locindex2extract(1,ind)
                npart = npartpergroup(groupid,ifile)
#ifdef DEBUG
                If(locindex2extract(2,ind) /= ifile) Then
                   Print *,'Pb: trying to read a group in file ',ifile, 
                   & ' when it should be found in file ',locindex2extract(2,ind)
                End If
#endif
                If( npart/=0 ) Then ! there is at least 1 part. to read in this group
                   Write(charigroup(1:8),'(I8.8)') groupid
                   grname='group'//charigroup
                   Call hdf5_open_group(dataid,grname,gr_id)
                
                   !read datasets
                   dname='position_part'
                   Call hdf5_read_data(gr_id, dname, 3, npart, pos(:,indpart:indpart+npart-1))

                   dname='velocity_part'
                   Call hdf5_read_data(gr_id, dname, 3, npart, vel(:,indpart:indpart+npart-1))

                   dname='identity_part'
                   Call hdf5_read_data(gr_id, dname, npart, id(indpart:indpart+npart-1))

                   If(potential) Then
                      dname='potential_part'
                      Call hdf5_read_data(gr_id, dname, npart, pot(indpart:indpart+npart-1))
                   End If

                   If(force) Then
                      dname='gravitational_field_part'
                      Call hdf5_read_data(gr_id, dname, 3, npart, for(:,indpart:indpart+npart-1))
                   End If

                   Call hdf5_close_group(gr_id)
                End If

                indpart = indpart+npart
             End Do
             Call hdf5_close_group(dataid)
             Call hdf5_close_file(file_id)
          End If
       End If
    End Do
    

    ! Don't forget to deallocate pos, vel, id and pot !
    Deallocate(boundaries, npartpergroup, globindexcube, coord2file)
    Deallocate(locindex2extract, ngroupperfile)
    
    ! If you use hdf5 later in your code, you can comment this
    Call hdf5_finalize() 

  End Subroutine extract_cuboid


  !=======================================================================
  !> Quick sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Recursive Subroutine quicksort(p,r,tab)
    
    Implicit None

    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted

    ! Input parameters modified in subroutine
    Integer(kind=4),Intent(inout),dimension(2,*)   :: tab   !< array

    ! Local parameters
    Integer(kind=4) :: q

#ifdef DEBUG2
    Print *,'quicksort begins'
#endif

    If(p<r) Then
       
       Call partition(p,r,q,tab)
       Call quicksort(p,q-1,tab)
       Call quicksort(q+1,r,tab)

    End If

#ifdef DEBUG2
    Print *,'quicksort ends'
#endif

  End Subroutine quicksort


  !=======================================================================

  Subroutine partition(p,r,q,tab)

    Implicit None

    ! Input parameters
    Integer(kind=4),Intent(in)  :: p, r  ! First and last index of the tables to sort
    ! Input parameters modified in subroutine
    Integer(kind=4),Intent(inout),dimension(2,*)   :: tab
    ! Output parameters
    Integer(kind=4), Intent(out) :: q
    ! Local parameters
    Integer(kind=4), dimension(2) :: tmpi
    Integer(kind=4) :: tref
    Integer(kind=4) :: i, j

#ifdef DEBUG2
    Print *,'partition begins'
#endif

    tref  = tab(2,r)

    i = p-1
    Do j = p, r-1
       If(tab(2,j) <= tref) Then
          i = i+1
          tmpi = tab(:,i)
          tab(:,i) = tab(:,j)
          tab(:,j) = tmpi
       End If
    End Do

    tmpi = tab(:,i+1)
    tab(:,i+1) = tab(:,r)
    tab(:,r) = tmpi

    q = i+1

#ifdef DEBUG2
    Print *,'partition ends'
#endif

  End Subroutine partition



End Module modextract