program sphericaloverdensity
  !---------------------------------------------------------------------
  ! Ce programme calcule la fonction de masse avec
  ! l'algorithme SPHERICAL OVERDENSITY.
  ! Ce programme doit lire en input les fichiers suivants:
  !          - un fichier resultat RAMSES:    part_00600.out
  ! Il genere en output les fichiers suivants:
  !          - un fichier fonction de masse:  mf_sod_00600.dat
  !          - un fichier structure:          st_sod_00600.dat
  ! E. Audit, R. Teyssier
  ! Meudon, le 30/05/96.
  !---------------------------------------------------------------------
  ! Version F90 par R. Teyssier le 28/02/00   
  ! f90 sod.f90 -o sod
  ! Memory: 44*npart in bytes
  !         npart=256**3   707 Mo
  !         npart=512**3   5.6 Go
  !---------------------------------------------------------------------
  ! Ajout du calcul de la masse pour un deuxieme seuil de sur densite
  ! Troisieme fichier d'output contenant cette deuxieme masse:
  !          - mass_sod_'num'.dat
  ! I. Balmes, 02/2012
  !---------------------------------------------------------------------
  ! Modification pour calculer le profil à partir des halos sod.
  ! Fichier d'output contenant pour chaque halo 30 bins de radius et 
  ! le nombre de particules correspondant.
  !          - halos_rad_part_'num'.dat
  ! I. Balmes, 01/2013
  !---------------------------------------------------------------------
  ! Ajout du module pour lire les fichiers fof, et transformation du code
  ! pour ne lire qu'un fichier fof (zone utile) 
  ! et ses voisins (zone buffer). 
  ! Le code recherche les halos dont le centre est dans la zone utile.
  ! version mpi, chaque processus lit plusieurs fichiers fof successivement.
  ! I. Balmes, 01/2014
  !---------------------------------------------------------------------
  ! Les outputs sont maintenant au format fof:
  ! fichier masst : contient pour chaque halo
  ! nombre de part, position du centre, numero du halo
  ! fichier strct : contient
  ! nb de halo
  ! puis pour chaque halos
  ! nb de part, positions de chaque part, vitesse de chaque part, 
  ! identifiant de chaque part
  ! compile with : mpif90 -O3 -g -traceback -limf  -convert big_endian sod_profile_rdfof_parallel.f90 -o sod_parallel
  ! I. Balmes, 02/2014
  !---------------------------------------------------------------------

  use mpi
  implicit none
  integer::ndim,npart,ngrid,n,i,icpu,ipos,ipart,nparttmp
  integer::ny,nz,ncpu,npart_new,nzone,nparttotal
  integer::ncpu2,npart2,ndim2
  integer::idm
  
  real(kind=8)::r,scale,Mmin,mtot,rhomoyen
  real,dimension(:,:),allocatable::xtmp,x,vel,vtmp,v
  real,dimension(:),allocatable::mp,mptmp
  real,dimension(:),allocatable::ap
  real(kind=8),dimension(:),allocatable::tmp
  real,dimension(:),allocatable::tmpmp
  real,dimension(:,:),allocatable::tmpx

  integer,dimension(:),allocatable::inttmp,idtmp
  integer,dimension(:),allocatable::isort,tmpid,idpreal8
  integer,dimension(:),allocatable::mask
  real(kind=8),dimension(:),allocatable::xsort

  real(kind=8)::seuil=600.,seuil2=800.
  integer::nmin=100,nx=0,ncut=4

  character*5::nchar,ncharcpu
  character*200::nomfich,root,repository
  logical::ok
  
  integer,parameter :: IRandNumSize = 4
  integer,dimension(IRandNumSize) :: localseed
  real(kind=8)::mstar_tot,mstar_lost
  integer::nsink,nstar_tot

  integer  ::nslice
  integer  ::j,nloc,proc_ID,k
  integer(8)::ll,kk!x,npart
  character::idp*5,filename*250

  real(kind=4), dimension(3)   :: mintmp,maxtmp
  real(kind=4), allocatable    :: cdm(:,:)
  integer(kind=8), allocatable :: id(:)
  integer(kind=8)              :: nthalo

  real(kind=4)::radius
  real(kind=4)::center(3)

  integer :: ixzone,iyzone,izzone,nxzone,nxbuffer
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  integer :: nxcpu,nycpu,nzcpu,itmp,jtmp,ktmp
  integer,dimension(27) :: list_cpu
  real :: xmin,xmax,ymin,ymax,zmin,zmax
  real :: xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone
  real :: xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter
  real :: xc,yc,zc,rxsmall,rysmall,rzsmall,rxlarge,rylarge,rzlarge
  real :: dx,dy,dz

  ! IO variables
  integer   :: uout,tot
  character :: conetype*1
  logical   :: found,writings

  integer :: procNB,procID,mpierr
  integer :: mynbfile,nmod,firstp,lastp

  call MPI_INIT(mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)

!  integer,parameter :: IRandNumSize = 4
!  integer,dimension(IRandNumSize) :: localseed
!  real(kind=8)::mstar_tot,mstar_lost
!  integer::nsink,nstar_tot
  nsink=0
  nstar_tot=0 

  !!!!!!!ATTENTION IL FAUDRA DISTINGUER ETOILES ET PARTICULES QUAND RUN HYDRO!!!!
  call read_params

  ndim=3
  nparttmp=0

 ! Count the number of cube files present
  ncpu=0
  found=.true.
  do while(found.eqv..true.)
     write(idp,'(I5.5)') ncpu
     filename=trim(root)//'_cube_'//idp
!     write(*,*) filename
     inquire(file=filename, exist=found)
     ncpu=ncpu+1   
  enddo

  ncpu=ncpu-1
  write(*,*) ncpu

  ! Count the number of particles
  nthalo=0
  do kk=0,ncpu-1
!  do kk=7,nslice-1
     write(idp,'(I5.5)') kk
     filename=trim(root)//'_cube_'//idp
     open(40,file=filename, status='Old', form='Unformatted')
     read(40) nloc
     close(40)
     nthalo=nloc+nthalo
  enddo
  nparttotal=nthalo

  !PARALLELIZED LOOP: firstp=local start, lastp=local end
  mynbfile = ncpu/procNB
  nmod=mod(ncpu,procNB)
  if(procID.le.nmod-1)then
     mynbfile=mynbfile+1
     firstp=procID*mynbfile
     lastp=(procID+1)*mynbfile-1
  else
     firstp=procID*mynbfile+nmod
     lastp=(procID+1)*mynbfile+nmod-1
  endif

! number of the fof file to read 
  do nzone=firstp,lastp
   ! corresponding indexes
     nxcpu=int(real(ncpu)**(1./3.))
     nycpu=int(real(ncpu)**(1./3.))
     nzcpu=int(real(ncpu)**(1./3.)) 

     ixzone=int(nzone/(real(nycpu*nzcpu)))
     iyzone=int((nzone-ixzone*nycpu*nzcpu)/real(nzcpu))
     izzone=nzone-ixzone*nycpu*nzcpu-iyzone*nzcpu

     nxzone=nx/nxcpu

     ixmin=ixzone-1
     ixmax=ixzone+1
     iymin=iyzone-1
     iymax=iyzone+1
     izmin=izzone-1
     izmax=izzone+1

     write(*,*) 'From', procID, 'Reading fof file: ',nzone
     write(*,*) 'From', procID, 'indexes ',ixzone,iyzone,izzone
  
     kk=0
     do i=ixmin,ixmax 
        do j=iymin,iymax
           do k=izmin,izmax
              itmp=i
              jtmp=j
              ktmp=k
              if (itmp.lt.0) itmp=nxcpu-1   
              if (jtmp.lt.0) jtmp=nycpu-1
              if (ktmp.lt.0) ktmp=nzcpu-1
              list_cpu(kk)=mod(ktmp,nzcpu)+mod(jtmp,nycpu)*nzcpu+mod(itmp,nxcpu)*nzcpu*nycpu  !list_cpu(0) !!!RY
              kk=kk+1
           enddo
        enddo
     enddo

     ! Count the number of particles
     nthalo=0
     do kk=0,26
        write(idp,'(I5.5)') list_cpu(kk)
        filename=trim(root)//'_cube_'//idp
        open(40,file=filename, status='Old', form='Unformatted')
        read(40) nloc
        close(40)
        nthalo=nloc+nthalo
     enddo
     nparttmp=nthalo

     write(*,*)'From', procID, 'Found ',nparttmp,' particles (including eventual stars)'
     write(*,*)'From', procID, 'Reading positions and masses...'
     allocate(xtmp(1:nparttmp,1:ndim),vtmp(1:nparttmp,1:ndim))
     allocate(mptmp(1:nparttmp),idtmp(1:nparttmp))
     if(nstar_tot > 0)allocate(ap(1:nparttmp))

     write(*,*)'From', procID, 'Allocating done'

     uout=20
     ll=1

     npart=0
     !the buffer zone is half the length of the working zone
     nxbuffer=nxzone/2. !nxbuffer is real RY !!!!
     !computation for the limits in xmin,xmax of the entire zone studied
     xmin=(ixzone-real(nxbuffer)/nxzone)/real(nxcpu)
     xmax=(ixzone+1+real(nxbuffer)/nxzone)/real(nxcpu)
     ymin=(iyzone-real(nxbuffer)/nxzone)/real(nycpu)
     ymax=(iyzone+1+real(nxbuffer)/nxzone)/real(nycpu)
     zmin=(izzone-real(nxbuffer)/nxzone)/real(nzcpu)
     zmax=(izzone+1+real(nxbuffer)/nxzone)/real(nzcpu)

     xmincenter=(ixzone-real(nxbuffer)/(2*nxzone))/real(nxcpu)
     xmaxcenter=(ixzone+1+real(nxbuffer)/(2*nxzone))/real(nxcpu)
     ymincenter=(iyzone-real(nxbuffer)/(2*nxzone))/real(nycpu)
     ymaxcenter=(iyzone+1+real(nxbuffer)/(2*nxzone))/real(nycpu)
     zmincenter=(izzone-real(nxbuffer)/(2*nxzone))/real(nzcpu)
     zmaxcenter=(izzone+1+real(nxbuffer)/(2*nxzone))/real(nzcpu)

     !Boundary conditions: these quantities can be below 0 or above 1

     write(*,*) 'From', procID, 'Buffer zone:', xmin,xmax,ymin,ymax,zmin,zmax

     !limits in xmin,xmax without the buffer zone
     xminzone=ixzone/real(nxcpu)
     xmaxzone=(ixzone+1)/real(nxcpu)
     yminzone=iyzone/real(nycpu)
     ymaxzone=(iyzone+1)/real(nycpu)
     zminzone=izzone/real(nzcpu)
     zmaxzone=(izzone+1)/real(nzcpu)

     !center of the zone
     xc=(xminzone+xmaxzone)/2.
     yc=(yminzone+ymaxzone)/2.
     zc=(zminzone+zmaxzone)/2.

     !distances between the center and the limits
     rxsmall=(0.5+real(nxbuffer)/nxzone)/real(nxcpu)
     rysmall=(0.5+real(nxbuffer)/nxzone)/real(nycpu)
     rzsmall=(0.5+real(nxbuffer)/nxzone)/real(nzcpu)

     rxlarge=1-rxsmall
     rylarge=1-rysmall
     rzlarge=1-rzsmall

     write(*,*) 'From', procID, 'Working zone:', xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone

     loadfile: DO kk=0,26

        write(idp,'(I5.5)') list_cpu(kk)
        filename=trim(root)//'_cube_'//idp
        !     write(*,*) filename

        open(uout,file=filename, status='Old', form='Unformatted')
        read(uout) nloc
        read(uout) proc_ID
        ! this is xmin,xmax for the cube
        read(uout)mintmp(1),maxtmp(1),mintmp(2),maxtmp(2),mintmp(3),maxtmp(3)

        allocate(cdm(3,nloc),v(3,nloc),id(nloc))
        read(uout) ((cdm(j,i),j=1,3),i=1,nloc)
        read(uout) ((v(j,i),j=1,3),i=1,nloc)
        read(uout) (id(i),i=1,nloc)

        do i=1,nloc
           do j=1,3
              xtmp(ll,j)=cdm(j,i)
              vtmp(ll,j)=v(j,i)
           enddo
           !dealing with boundary conditions
           dx=abs(xtmp(ll,1)-xc)
           dy=abs(xtmp(ll,2)-yc)
           dz=abs(xtmp(ll,3)-zc)
           !this line counts the number of particles that will be kept
           if (dx.le.rxsmall.or.dx.ge.rxlarge.and.&
                dy.le.rysmall.or.dy.ge.rylarge.and.&
                dz.le.rzsmall.or.dz.ge.rzlarge) npart=npart+1
           mptmp(ll)=1.
           idtmp(ll)=id(i)
           ll=ll+1
           !        write(*,*) ll
           !if (kk.eq.1025) then 
           !   write(*,*) i,'/',nloc
           !end if
        enddo
        !     write(*,*) 'cdm, vel and idtot have values'

        deallocate(cdm,v,id)
        close(uout)

     ENDDO loadfile

     write(*,*) 'From', procID, 'Found ', npart,' particles in the buffer zone'

     !too many particles have been read, lots need to be thrown away
     allocate(mask(1:npart))
     allocate(x(1:npart,1:ndim),vel(1:npart,1:ndim))
     allocate(idpreal8(1:npart),mp(1:npart))
     allocate(xsort(1:npart),isort(1:npart))
     ll=1
     do i=1,nparttmp                                                
        dx=abs(xtmp(i,1)-xc)
        dy=abs(xtmp(i,2)-yc)
        dz=abs(xtmp(i,3)-zc)
        if (dx.le.rxsmall.or.dx.ge.rxlarge.and.&
             dy.le.rysmall.or.dy.ge.rylarge.and.&
             dz.le.rzsmall.or.dz.ge.rzlarge) then
           mask(ll)=i
           ll=ll+1
        endif
     enddo
     
     x=xtmp(mask,:)    
     vel=vtmp(mask,:)
     idpreal8=idtmp(mask)
     mp=mptmp(mask)

     deallocate(xtmp,vtmp,idtmp,mptmp,mask)
     if(nstar_tot > 0)deallocate(ap)

     write(*,*) 'From', procID, 'particles out of the buffer zone have been thrown away'

     !we actually donot care about velocities nor ids (?!)
     !  deallocate(vel,id)

     !on suppose qu'il n'y a pas d'étoiles

     !nx must be initialized in the arguments. indicate that it is not optional anymore!

     ny=nx
     nz=nx
     rhomoyen=1.0
     write(*,*)'Working mesh=',nx,ny,nz
     write(*,*)'From', procID, 'Found rho_bar=',rhomoyen,' -> rescaling to 1.0'
     mp=1.0 !useless RY !!!

     scale=nx
     x=x*scale
     xminzone=xminzone*scale
     xmaxzone=xmaxzone*scale
     yminzone=yminzone*scale
     ymaxzone=ymaxzone*scale
     zminzone=zminzone*scale
     zmaxzone=zmaxzone*scale
     !this is just a factor 1.
     Mmin=dble(nmin)!/dble(npart)*scale**3


     ngrid=nx*ny*nz

!     write(*,*) minval(x),maxval(x)

     !------------------------------------------
     ! Tri des particules selon leur densite NGP
     !------------------------------------------
     write(*,*) 'From', procID, 'Sorting particles by NGP density'
     call tri_ngp(x,mp,xsort,isort,npart,nx,ny,nz,npart_new)

!     write(*,*) minval(x),maxval(x)

     !---------------------------------------------------------
     ! Calcul de la fonction de masse par SPHERICAL OVERDENSITY
     !---------------------------------------------------------
     write(*,*) 'From', procID, 'Computing halo mass function'
     write(*,*) 'using a spherical overdensity delta_bar=',seuil
     write(*,*) 'using a second overdensity delta_bar=',seuil2
     call spherover(x,vel,mp,idpreal8,xsort,isort,nchar,npart,nx,ny,nz,&
          npart_new,ncut,seuil,seuil2,Mmin,scale,nparttotal,nzone,&
          xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone,&
          xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter)
     write(*,*) 'From', procID, ' done ',nzone
     deallocate(x,vel,idpreal8,mp,xsort,isort)
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  if(procID==0)write(*,*) 'Run completed!'
  call MPI_FINALIZE(mpierr)

  stop

contains
  subroutine read_params

      implicit none

      integer       :: i,n
      integer       :: iargc
      character(len=4)   :: opt
      character(len=128) :: arg
      LOGICAL       :: bad, ok
      
      n = iargc()
      if (n < 2) then
         print *, 'usage: sod  [-inp input_dir]'
         print *, '            [-nx  nx_grid]  (optional)'
         print *, '            [-min np_min]   (optional)'
         print *, '            [-sel ncell]    (optional)'
         print *, '            [-del del_bar]  (optional)'
         print *, '            [-dl2 del_bar2](optional)'
         print *, '            [-hlp]          (optional)'
         call MPI_FINALIZE(mpierr)
         stop
      end if

      do i = 1,n,2
         call getarg(i,opt)
         if (opt == '-hlp') then
            print 1, repository,nx,seuil,mmin
            call MPI_FINALIZE(mpierr)
            stop
         end if
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            call MPI_FINALIZE(mpierr)
            stop 2
         end if
         call getarg(i+1,arg)
         select case (opt)
         case ('-inp')
            repository = trim(arg)
         case ('-nx')
            read (arg,*) nx
         case('-min')
            read (arg,*) nmin
         case('-sel')
            read (arg,*) ncut
         case ('-del') 
            read (arg,*) seuil
         case ('-dl2')
            read (arg,*) seuil2
         case ('-nch')
            read (arg,*) nchar
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      root=repository

1     format(/,&
           & " -inp [] ",A,/,&
           & " -nx  [] ",I6,/,&
           & " -min [] ",I6,/,&
           & " -sel [] ",I6,/,&
           & " -del [] ",e10.2,/)

      return

    end subroutine read_params

end program sphericaloverdensity

!------------------------------------------------------------
! TRI DES PARTICULES SELON LA DENSITE NGP SUR UNE GRILLE FINE
!------------------------------------------------------------
subroutine tri_ngp(x,mp,pdens,isort,npart,nx,ny,nz,npart_new)
  use mpi
  implicit none
  integer::mpierr,procID,procNB
  integer::npart,nx,ny,nz,npart_new
  integer,dimension(1:npart)::isort
  real,dimension(1:npart,1:3)::x
  real,dimension(1:npart)::mp
  real(kind=8),dimension(1:npart)::pdens

  integer(kind=8),dimension(:),allocatable::indx
  integer::ngxngy,ngx,ngy,ngz
  integer(kind=8)::nmax,i1,i2,i3,ind
  integer::i,j,ifail
  integer::imin,imax,ntot
  real(kind=8)::amass

  ngx=4*nx
  ngy=4*ny
  ngz=4*nz

  ngxngy=ngx*ngy

  write(*,*) 'npart here is ',npart
   
  allocate(indx(1:npart))

  ! Remplissage du tableau des indices
  !-----------------------------------
  do i = 1,npart
     i1 = int(x(i,1)*dble(ngx)/dble(nx))
     i2 = int(x(i,2)*dble(ngy)/dble(ny))
     i3 = int(x(i,3)*dble(ngz)/dble(nz))
     indx (i) = 1+i1+i2*ngx+i3*ngxngy
  end do
  
  ! Tri des particules selon les indices
  !-------------------------------------
  ifail=0
  pdens=indx
  call quick_sort(pdens,isort,npart)
  indx=pdens

  ! Calcul de l'inverse de la densite NGP pour chaque particule
  !------------------------------------------------------------
  imin = 1
  do while(imin .le. npart)
     ind  = indx(imin)
     i    = imin
     do while(indx(i).eq.ind.and.i.lt.npart)
        i = i + 1
     end do
     if(i.eq.npart)i=npart+1
     imax = i - 1
     ntot = i - imin
     amass= 0.0
     do i = imin ,imax
       amass=amass+mp(isort(i))
     end do
     do i = imin ,imax
        pdens(isort(i)) = 1.0d0/amass
     end do
     imin = imax + 1
  end do
  print *, 1./minval(pdens)
  print *, 1./maxval(pdens)

  ! Tri des particules selon la densite
  !------------------------------------
  call quick_sort(pdens,isort,npart)

  i=1
  do while (pdens(i) < 0.5)
     i=i+1
  end do
  npart_new=i
  write(*,*)'npart_active=',npart_new
  deallocate(indx)
  return

end subroutine tri_ngp

!---------------------------------------------------------------------------
subroutine spherover(x,vel,mp,idpreal8,xsort,isort,nchar,npart,nx,ny,nz,&
                     npart_new,ncut,seuil,seuil2,Mmin,scale,nparttotal,nzone,&
                     xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone,&
                     xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter)
!---------------------------------------------------------------------------
  use mpi
  implicit none
  integer::mpierr,procID,procNB
  integer,parameter:: nbins=60
  character*5::nchar
  integer::npart,nx,ny,nz,npart_new,nparttotal
  real(kind=8)::seuil,Mmin,scale,seuil2
  integer::ncut,ncut2          ! Amplitude de selection des particules
  real,dimension(1:npart,1:3)::x,vel
  real,dimension(1:npart)::mp
  real(kind=8),dimension(1:npart)::xsort
  integer,dimension(1:npart)::isort,idpreal8
  real::xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone  
  real::xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter  

  integer::np,nzone
  real(kind=8)::twopi=6.283185307179586
  real(kind=8)::dtmax=0.01 ! dr/rayon pour la convergence du barycentre
  integer::Nsel=4000000    ! Nombre maximum de particules selectionnees

  integer,dimension(:),allocatable::structure,allpart
  integer,dimension(:),allocatable::strct,masses
  integer,dimension(:),allocatable::indx_ngp,isort_ngp
  integer,dimension(:),allocatable::first_ngp,num_ngp,idp

  integer,dimension(:),allocatable::indx_sel,isort_sel,irank_sel
  real,dimension(:),allocatable::m_sel
  real,dimension(:,:),allocatable::x_sel,x_sel2,v_sel
  real(kind=8),dimension(:),allocatable::distance
  real(kind=8)::rayon_min,distance_new

  integer::i,j,k,l,npart_sel,ipart,npart_sel2,n_max,n_test
  integer::is,iamas,namas,gap,masse,ip,ifail,masse2
  integer::ic,jc,kc,ii,jj,kk,ind,indc
  integer::i1,i2,i3,imin,imax,ntot
  integer::niter,itest
  real(kind=8)::xc,yc,zc,xx,yy,zz,xc0,yc0,zc0,dc0
  real(kind=8)::xcnew,ycnew,zcnew,xxnew,yynew,zznew,dxnew,dynew,dznew
  real(kind=8)::xb,yb,zb,dt
  real(kind=8)::taillex,tailley,taillez,amass
  real(kind=8)::dx,dy,dz,d,d2,rayon,rayon2,volume,overdens,overdens2
  real(kind=8)::dmin,x_tmp,y_tmp,z_tmp,dis_tmp,qqq
  real(kind=8)::Mtotsel

  character*5::nzonechar
  character*80::nomf,nomf2,filestrct,fileamas,fileprofile,filestrctfof,filemasst
  integer::count_no_halo
  real::frac=1 !!! Control when to stop the loop on particles. Added by Yann ! a diminuer ?

  real(kind=8),dimension(1:nbins)::rad_bins
  integer,dimension(1:nbins)::mass_bins

!  write(*,*) minval(x),maxval(x)

  allocate(structure(1:npart),strct(1:nparttotal),indx_ngp(1:npart),isort_ngp(1:npart))
  allocate(masses(1:npart),allpart(1:npart))
! masses contains the masses of all the halos found. 
! allpart contains the index of the parts of each amas (in x, vel and idpreal8 arrays) (from 1 to masses(1), ids of part in first amas, etc)

  structure = 0
  taillex   = dble(nx)
  tailley   = dble(ny)
  taillez   = dble(nz)
  rayon_min = sqrt(3.)*scale/(2.*nx)
  write(nzonechar,'(I5.5)') nzone
  nomf      = 'sod_'//TRIM(nchar)//'_'//TRIM(nzonechar)//'.dat'
  nomf2     = 'mass_sod_'//TRIM(nchar)//'_'//TRIM(nzonechar)//'.dat'
  fileamas  = 'm_sod_'//trim(nchar)//'_'//TRIM(nzonechar)//'.dat'
  fileprofile = 'halos_rad_part_'//trim(nchar)//'_'//TRIM(nzonechar)//'.dat'
  filestrctfof = 'sod_strct_'//trim(nchar)//'_'//trim(nzonechar)//'.dat'
  filemasst = 'sod_masst_'//trim(nchar)//'_'//trim(nzonechar)//'.dat'

  allocate(num_ngp  (1:nx*ny*nz))
  allocate(first_ngp(1:nx*ny*nz))

  allocate(indx_sel (Nsel))
  allocate(isort_sel(Nsel))
  allocate(distance (Nsel))
  allocate(x_sel    (Nsel,1:3))
  allocate(v_sel    (Nsel,1:3))
  allocate(x_sel2    (Nsel,1:3))
  allocate(m_sel    (Nsel))
  allocate(idp      (Nsel))
  ncut2=1

  ! Calcul du pointeur cellules vers particules
  !--------------------------------------------
  write(*,*)'Building particle linked list...'
  do i  = 1,npart
     i1 = int(x(i,1))
     i2 = int(x(i,2))
     i3 = int(x(i,3))
     indx_ngp(i) = 1+i1+i2*nx+i3*nx*ny
  end do

  xsort=indx_ngp
  call quick_sort(xsort,isort_ngp,npart)
  indx_ngp=xsort

! as far as I understand
! first_ngp contains the index of the first particule in cell i
! num_ngp contains the number of part in cell i
! indx_ngp contains the number of the cell in which is the part
  first_ngp(:)=0
  num_ngp(:)=0
  imin = 1
  do while(imin.le.npart)
     ind  = indx_ngp(imin)
     i    = imin
     first_ngp(ind) = imin
     do while(indx_ngp(i).eq.ind.and.i.lt.npart)
        i = i + 1
     end do
     if(i.eq.npart)i=npart+1
     imax = i - 1
     ntot = i - imin
     num_ngp(ind) = ntot
     imin = imax + 1
  end do

  write(*,*) '  halo#    npart   npart2   mass       x          y          z         it  radius     del   del2'
  open(10,file=nomf,form='formatted')
  open(20,file=nomf2,form='formatted')
  open(22,file=fileamas,form='formatted')
  open(25,file=fileprofile,form='formatted')
  open(26,file=filestrctfof,form='unformatted')
  open(27,file=filemasst,form='unformatted')
  namas = 0
  gap=0
  xsort = 1.0
  count_no_halo=0

  ipart=1
  !do ipart = 1,npart_new
  do while( (ipart.le.npart_new).and.(dble(count_no_halo)/dble(npart_new).le.frac))
  !do while( (ipart.le.npart_new))
     
     ! On choisit la particule la plus dense de la liste
     !-------------------------------------------------
     if(xsort(isort(ipart)).gt.0.)then
        
        ! On initialise le barycentre de l'amas avec 
        ! la position de la particule courante
        xc = x(isort(ipart),1)
        yc = x(isort(ipart),2)
        zc = x(isort(ipart),3)
! the halo identifier is the id of the barycentre
        iamas = idpreal8(isort(ipart))

        xc0=xc; yc0=yc; zc0=zc
        
        ! Preselection des particules les plus proches du barycentre
        !-----------------------------------------------------------
        is = 0
        ic = int(xc)
        jc = int(yc)
        kc = int(zc)
        indc=1+ic+jc*nx+kc*nx*ny
        mtotsel=0.
        do i=-ncut,ncut
           ii=i+ic
           if(ii.ge.nx)ii=ii-nx
           if(ii.lt.0) ii=ii+nx
           do j=-ncut,ncut
              jj=j+jc
              if(jj.ge.ny)jj=jj-ny
              if(jj.lt.0) jj=jj+ny
              do k=-ncut,ncut
                 kk=k+kc
                 if(kk.ge.nz)kk=kk-nz
                 if(kk.lt.0) kk=kk+nz
                 ind = 1+ii+jj*nx+kk*nx*ny
                 imin = first_ngp(ind)
                 imax = imin+num_ngp(ind)-1
                 do ip = imin,imax
                    if(xsort(isort_ngp(ip)).gt.0.)then
                       is = is + 1
                       if(is > Nsel)then
                          write(*,*)'WARNING: from ',nzonechar,' Increase Nsel in main program'
                          deallocate(structure,strct,indx_ngp,isort_ngp)
                          deallocate(masses,allpart)
                          deallocate(num_ngp,first_ngp,indx_sel)
                          deallocate(isort_sel,distance,x_sel,v_sel)
                          deallocate(x_sel2,m_sel,idp)
                          return
                       end if
                       indx_sel(is) =  isort_ngp(ip)
                       x_sel(is,1) = x(isort_ngp(ip),1)
                       x_sel(is,2) = x(isort_ngp(ip),2)
                       x_sel(is,3) = x(isort_ngp(ip),3)
                       v_sel(is,1) = vel(isort_ngp(ip),1)
                       v_sel(is,2) = vel(isort_ngp(ip),2)
                       v_sel(is,3) = vel(isort_ngp(ip),3)
                       m_sel(is) = mp (isort_ngp(ip))
                       idp(is) = idpreal8(isort_ngp(ip)) !index of the particule
                       mtotsel = mtotsel + m_sel(is)
                    end if
                 end do
              end do
           end do
        end do
        npart_sel = is
!        write(*,*)'npart_sel=',npart_sel
!        write(*,*)'Mtotsel=',Mtotsel,mmin

        !select particules that are in this cell plus 1 cell around for the barycentre loop
        is = 0
        do i=-ncut2,ncut2
           ii=i+ic
           if(ii.ge.nx)ii=ii-nx
           if(ii.lt.0) ii=ii+nx
           do j=-ncut2,ncut2
              jj=j+jc
              if(jj.ge.ny)jj=jj-ny
              if(jj.lt.0) jj=jj+ny
              do k=-ncut2,ncut2
                 kk=k+kc
                 if(kk.ge.nz)kk=kk-nz
                 if(kk.lt.0) kk=kk+nz
                 ind = 1+ii+jj*nx+kk*nx*ny
                 imin = first_ngp(ind)
                 imax = imin+num_ngp(ind)-1
                 do ip = imin,imax
                    if(xsort(isort_ngp(ip)).gt.0.)then
                       is = is + 1
                       x_sel2(is,1) = x(isort_ngp(ip),1)
                       x_sel2(is,2) = x(isort_ngp(ip),2)
                       x_sel2(is,3) = x(isort_ngp(ip),3)
                    end if
                 end do
              end do
           end do
        end do
        npart_sel2 = is

        if(Mtotsel .ge. Mmin) then

        dt    = 1.
        niter = 0
        

! debut boucle barycentre (?)
! 02/01/2013 : suppression de la boucle sur le barycentre.
! 14/01/2013 : boucle sur les particules voisines de celle de départ,
! pour déterminer celle qui a le plus de voisins.
!        do while(dt.gt.dtmax.and.niter.le.100)

           dx = abs(xc-xc0)
           dx = min(dx,taillex-dx)
           dy = abs(yc-yc0)
           dy = min(dy,tailley-dy)
           dz = abs(zc-zc0)
           dz = min(dz,taillez-dz)
           dc0 = dx*dx+dy*dy+dz*dz

           ! Calcul des distances des particules au barycentre
           !--------------------------------------------------
           do i = 1,npart_sel
              xx = x_sel(i,1)
              yy = x_sel(i,2)
              zz = x_sel(i,3)
              dx = abs(xx-xc)
              dx = min(dx,taillex-dx)
              dy = abs(yy-yc)
              dy = min(dy,tailley-dy)
              dz = abs(zz-zc)
              dz = min(dz,taillez-dz)
              distance(i) = dx*dx+dy*dy+dz*dz
           end do
          
           ! Tri des distances
           !---------------------------------------------------
           call quick_sort(distance,isort_sel,npart_sel)

           ! Boucle sur les particules les plus proches du barycentre
           ! pour rechercher la plus dense
           !---------------------------------------------------
           !loop on all particules in the same cell as the current barycentre
           n_max=0
           l=first_ngp(indc)
           !write(*,*) l, npart, nsel, indc, indx_ngp(l)
           !write(*,*) 'barycentre before loop :', xc,yc,zc
           do while(indx_ngp(l).eq.indc.and.l.lt.npart)
              !write(*,*) 'entered loop'
              !x(isort(ipart) ?
              xcnew=x(isort_ngp(l),1)
              ycnew=x(isort_ngp(l),2)
              zcnew=x(isort_ngp(l),3)
              !write(*,*) 'candidate :',xcnew/scale,ycnew/scale,zcnew/scale,taillex/scale,nx

              !compute the distances to the new barycentre for the preselected particules
              ! and count the number of particules in the sphere of radius rayon_min
              i=0
              do j=1,npart_sel2
                 xxnew=x_sel2(j,1)
                 yynew=x_sel2(j,2)
                 zznew=x_sel2(j,3)
                 dxnew=abs(xxnew-xcnew)
                 dxnew=min(dxnew,taillex-dxnew)
                 dynew=abs(yynew-ycnew)
                 dynew=min(dynew,tailley-dynew)
                 dznew=abs(zznew-zcnew)
                 dznew=min(dznew,taillez-dznew)
                 distance_new=dxnew*dxnew+dynew*dynew+dznew*dznew
                 !write(*,*) distance_new, rayon_min**2
                 if (distance_new.le.rayon_min**2) then
                    i=i+1
                 end if
              end do
              n_test=i
              !write(*,*) n_test, n_max
              !update the barycentre if there are more particules in the sphere
              if (n_test.ge.n_max) then
                 n_max=n_test
                 xc=xcnew
                 yc=ycnew
                 zc=zcnew
                 iamas=idpreal8(isort_ngp(l))
                 !write(*,*) 'new barycentre :', xc,yc,zc
              end if
              l=l+1
           end do


           !write(*,*) 'barycentre after loop :', xc,yc,zc
           !write(*,*) n_max,npart_sel2

           ! Calcul des distances des particules au barycentre
           !--------------------------------------------------
           do i = 1,npart_sel
              xx = x_sel(i,1)
              yy = x_sel(i,2)
              zz = x_sel(i,3)
              dx = abs(xx-xc)
              dx = min(dx,taillex-dx)
              dy = abs(yy-yc)
              dy = min(dy,tailley-dy)
              dz = abs(zz-zc)
              dz = min(dz,taillez-dz)
              distance(i) = dx*dx+dy*dy+dz*dz
           end do
          
           ! Tri des distances
           !---------------------------------------------------
           call quick_sort(distance,isort_sel,npart_sel)

           ! Calcul du rayon de Viriel et de la masse de l'amas
           !---------------------------------------------------
           i        = 1
           overdens = 2.0*seuil
           amass    = m_sel(isort_sel(i))
           rayon    = 0.0

           !sqrt or not sqrt ?
!           do while(i.lt.10.or.(overdens.gt.seuil.and.i.lt.npart_sel.and.rayon.lt.(real(ncut)-dc0)))
           do while(i.lt.10.or.(overdens.gt.seuil.and.i.lt.npart_sel.and.rayon.lt.(real(ncut)-sqrt(dc0))))
              i        = i+1
              rayon    = SQRT(distance(i))
              volume   = 2./3.*twopi*rayon**3
              amass    = amass+m_sel(isort_sel(i))
              overdens = amass/volume
!              write(*,*)i,rayon,amass,overdens
              ! Ajout du calcul de la masse de l'amas avec un deuxieme seuil de densite (seuil2>seuil)
              ! 02/2012, I.Balmes
              !---------------------------------------------------------------------------------------
              if (overdens.gt.seuil2) then 
                 masse2=i
                 if (rayon>1e-30) then
                    overdens2=overdens
                 else
                    overdens2=0
                 endif
              end if
           end do

           if(overdens.gt.seuil)then
              write(*,*)'WARNING: from ',nzonechar,' Increase selection radius (-sel option in command line)'
              write(*,*)'mass     =',i
              write(*,*)'nsel     =',npart_sel
              write(*,*)'overdens =',overdens
              write(*,*)'radius   =',rayon+dc0
              write(*,*)'nsel     =',ncut
              write(*,*)'niter    =',niter
              write(*,*)'mass2    =',masse2
              write(*,*)'overdens2=',overdens2
              deallocate(structure,strct,indx_ngp,isort_ngp)
              deallocate(masses,allpart)
              deallocate(num_ngp,first_ngp,indx_sel,isort_sel)
              deallocate(distance,x_sel,v_sel,x_sel2,m_sel,idp)
              return
           else
              rayon =SQRT(distance(i-1))
              volume=2./3.*twopi*rayon**3
              amass =amass-m_sel(isort_sel(i))
              if(volume > 0.)then
                 overdens=amass/volume
              else
                 overdens=0.
              end if
              masse=i-1
           end if
!           write(*,*)'iter=',niter,'mass=',masse

! 02/01/2013 : ajout d'une boucle pour le calcul du profil M(<r)
           do j=1,nbins
              !bins de rayon repartis logarithmiquement entre le rayon et le rayon*0.01
              rad_bins(j)=rayon*exp(log(0.01)+(j-1)*(log(1.)-log(0.01))/(nbins-1.))
              if (j.eq.1) then
                 k=0
              else 
                 k=mass_bins(j-1)
              end if
              do while(sqrt(distance(k)).lt.rad_bins(j))
                 mass_bins(j)=k
                 k=k+1
              end do
           end do      

!the halo is kept if its mass is above Mmin and its center is not in the buffer zone
        if(amass.ge.Mmin.and.xc.le.xmaxzone.and.xc.gt.xminzone&
                        .and.yc.le.ymaxzone.and.yc.gt.yminzone&
                        .and.zc.le.zmaxzone.and.zc.gt.zminzone)then
           namas=namas+1
           masses(namas)=amass
           gap=gap+amass
           write(*,997)iamas,masse,masse2,dble(amass)/scale**3,xc/scale,yc/scale,zc/scale &
                 & ,niter,rayon/scale,overdens,overdens2
           write(20,998)iamas,masse,masse2,overdens,overdens2
           write(22,996)iamas,masse
! ecriture du profil du halo M(<r)
           write(25,995)iamas,masse,xc/scale,yc/scale,zc/scale
           do i=1,nbins
              write(25,994)rad_bins(i)/scale,mass_bins(i)
           end do
           write(27)masse,xc/scale,yc/scale,zc/scale,iamas
           count_no_halo=0
        end if

        ! On retire de la liste des particules celles de l'amas
        !------------------------------------------------------
        ! version parallele
        ! on ne retire pas les particules de la liste si le 
        ! centre du halo est trop proche du bord de la zone buffer.
        !------------------------------------------------------
        do i = 1,masse
           if(amass.ge.Mmin.and.xc.le.xmaxzone.and.xc.gt.xminzone&
                           .and.yc.le.ymaxzone.and.yc.gt.yminzone&
                           .and.zc.le.zmaxzone.and.zc.gt.zminzone)then
              structure(indx_sel(isort_sel(i))) = iamas
              strct(idp(isort_sel(i)))=iamas
              allpart(i+gap-masse)=isort_sel(i)
              xsort(indx_sel(isort_sel(i))) = 0.
           else !if (xc.le.xmaxcenter.and.xc.gt.xmincenter&
               !.and.yc.le.ymaxcenter.and.yc.gt.ymincenter&
               !.and.zc.le.zmaxcenter.and.zc.gt.zmincenter)then
              structure(indx_sel(isort_sel(i)))=0
              strct(idp(isort_sel(i)))=0
              xsort(indx_sel(isort_sel(i))) = 0.
           endif
        enddo
!        write(*,*) 'Particules taken out'
	      
        ! Ecriture des caracteristiques de l'amas sur fichier
        !----------------------------------------------------
        if(amass.ge.Mmin.and.xc.le.xmaxzone.and.xc.gt.xminzone.and.&
                             yc.le.ymaxzone.and.yc.gt.yminzone.and.&
                             zc.le.zmaxzone.and.zc.gt.zminzone)then
           write(10,999)iamas,masse,amass/scale**3,xc/scale,yc/scale,zc/scale,rayon/scale
        end if
        
     end if
     
  end if
  count_no_halo=count_no_halo+1

  !if(mod(ipart,10000)==0)write(*,'(0PF5.1,"% complete")')100.*dble(ipart)/dble(npart_new)
  
  if(mod(ipart,npart_new/10)==0)then   
     write(*,'(0PF5.1,"% complete")')100.*dble(ipart)/dble(npart_new)  
  endif
  
  ipart=ipart+1
end do
   
  close(10)
  close(20)
  
  if(dble(count_no_halo)/dble(npart_new).gt.frac)then
     write(*,*)'Seems that there are no more halos... It s however not guaranteed! '
     write(*,*)'Increase the stopping criterion frac=',frac,' if you think there are still some more '
  endif

  nomf = 'st_sod_'//TRIM(nchar)//'.dat'
  open (10,file=nomf,form='unformatted')
  write(10)structure
  close(10)
  
  filestrct='strct_sod_'//trim(nchar)//'.dat'
  open(11,file=filestrct,form='unformatted')
  write(11) strct
  close(11)

  gap=0
  write(26) namas
  do i=1,namas 
     write(26) masses(i)
     do j=1,masses(i)
        write(26) x_sel(allpart(j),1)/scale,x_sel(allpart(j),2)/scale,x_sel(allpart(j),3)/scale
     enddo
     do j=1,masses(i)
        write(26) v_sel(allpart(j),1),v_sel(allpart(j),2),v_sel(allpart(j),3)
     enddo
     do j=1,masses(i)
        write(26) idp(allpart(j))
     enddo
     gap=gap+masses(i)
  enddo

  deallocate(structure,strct,indx_ngp,isort_ngp,num_ngp,first_ngp)
  deallocate(indx_sel,isort_sel,distance,x_sel,x_sel2,m_sel,idp)
  deallocate(v_sel,allpart,masses)

997 format (I12,1x,I8,1x,I8,4(1x,1pe10.3),1x,I2,1x,6(E10.3,1x))
998 format (I12,2(1x,I8),2(1x,E10.3))
999 format (I12,1x,I8,5(1x,1pe10.3))
996 format (I12,1x,I12)
995 format (I12,1x,I8,3(1x,1pe10.3))
994 format (E10.3,1x,I8)
	
  return
end subroutine spherover

subroutine title(n,nchar)
  use mpi
  integer::mpierr,procID,procNB
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5


  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif
end subroutine title

SUBROUTINE quick_sort(list, order, n)
  
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  USE mpi
  IMPLICIT NONE
  INTEGER :: n,mpierr,procID,procNB
  REAL*8, DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
  
  ! Local variable
  INTEGER :: i
  
  DO i = 1, n
     order(i) = i
  END DO
  
  CALL quick_sort_1(1, n)
  
CONTAINS
  
  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
    use mpi
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL*8              :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6
    
    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)
       
    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       
       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO
          
          
          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       
       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1
  
  
  SUBROUTINE interchange_sort(left_end, right_end)
    use mpi
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL*8              :: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO
    
  END SUBROUTINE interchange_sort
  
END SUBROUTINE quick_sort
