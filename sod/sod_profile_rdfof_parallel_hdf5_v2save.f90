module modvar
#ifdef LONGINT
  integer,parameter::idprec=8,idprecread=8
#else
  integer,parameter::idprec=4,idprecread=4
#endif

  integer(idprec),dimension(:),allocatable::idf !id
  real,dimension(:,:),allocatable::posf,velf   !position velocity
  real,dimension(:),allocatable::mpf   !mass

  real(8),dimension(:),allocatable::xsort !density
  integer,dimension(:),allocatable::isort !sort index
  integer::IOGROUPSIZE=0,dummy_io,tag=2700,info  !IO tickets parameters

end module modvar

module modcic
  real(8),dimension(:,:,:),allocatable::cube
  real(8),dimension(:),allocatable::rho
end module modcic

program sphericaloverdensity
  use modvar
  use modcic
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
  ! compile with : mpif90 -O3 -g -traceback -limf  -fp-model strict sod_profile_rdfof_parallel.f90 -o sod_profile_rdfof_parallel
  ! or for reading  hdf5 cube
  !  mpif90 -O3 -g -traceback -limf  -fp-model strict -cpp -DREADHDF5 -DLONGINT modhdf5.f90 modvariables.f90 modextract.f90 sod_profile_rdfof_parallel_hdf5.f90 -o sod_profile_rdfof_parallel_hdf5
  ! I. Balmes, 02/2014
  !Run with 
  !mpirun -np 64 ~/proj/quint/code/sod/sod_profile_rdfof_parallel 
  ! -inp ../../fof/output_00014/fof_boxlen162_n128_lcdmw5  -nx 128 -min 100 -sel 4 -del 200. -dl2 600. -nch 00014 -nxb 16 -nse 200000
  ! Heavy tests, debugging, cleaning + parallel CIC, Y.Rasera 09/2015
  ! remark: 03 only not 100% safe with openmpi+ifort need to put -fp-model strict
  !
  ! TO DO: adapt to hydro run + improve I/O by reading HDF5 with many subcubes inside each cube or
  ! one cpu read one cube and send to the other (like pFoF) + split main in several routines+translate above
  ! + do not restrict to one cube and 26 neighbours+remove unused variable + do not remove part in halo near edge
  ! of buffer + make profile written with io ticket +double check conversion integer real8 ok + test larger simu + link with pfof and compute halo prop
  !+ test version ultra original sod + test profile + selection moins large rayon_min + test occupation memoire etc... + integer8
  !---------------------------------------------------------------------

  use mpi
#ifdef READHDF5
  use modhdf5
  use modwritehalo
  use modextract
  Use modmpicom,only:procNB,procID

#endif

  implicit none

  integer::ndim,npart,i,ipos,ipart
  integer::ny,nz,ncpu,npart_new,nzone
  integer::ncpu2,npart2,ndim2,npartini
  real(8)::r,scale,Mmin,mtot

  real,dimension(:,:),allocatable::v,cdm
  integer(idprecread),dimension(:), allocatable :: idtmp

  integer::nx=0
  character*5::nchar=''
  
  real(8)::seuil=200.,seuil2=800.
  integer::nmin=100

  integer::ncut=4,nxbuffer=0,Nsel=4000000
  integer::nparthalo_upperbound_estimate_default=0 !if >0: try better estimate of
  !default values ncut, nxbuffer,nsel

  integer::nparthalo_upperbound_estimate=0
  
  integer::nthalo,nparttotal
  
  character*5::ncharcpu
  character*200::nomfich,root,repository
  
  integer  ::nslice
  integer  ::j,nloc,proc_ID,k,ifillpart
  integer(8)::kk
  character::idpchar*5,filenametmp*250


  real(kind=4), dimension(3)   :: mintmp,maxtmp

  integer :: ixzone,iyzone,izzone,nxzone
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  integer :: nxcpu,nycpu,nzcpu,itmp,jtmp,ktmp
  integer,dimension(0:26) :: list_cpu 
  real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(8) :: xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone
  real(8) :: xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter
  real(8) :: xc,yc,zc,rxsmall,rysmall,rzsmall,rxlarge,rylarge,rzlarge
  real(8) :: dx,dy,dz

  ! IO variables
  integer   :: uinout=20
  logical   :: found

  !integer :: procNB,procID
  integer::mpierr
  integer :: mynbfile,nmod,firstp,lastp

    !CIC variables
  logical:: docic=.true.
  real(8)   ::xxmin,xxmax,yymin,yymax,zzmin,zzmax 
  real(8)   ::xmincic,xmaxcic,ymincic,ymaxcic,zmincic,zmaxcic 
  integer::local_nx,local_ny,local_nz          
  integer::ixzonecic,iyzonecic,izzonecic       
  integer::nxzonecic,nyzonecic,nzzonecic       
  logical::periodic                            
  integer::nxcic,nycic,nzcic
  integer::nxngp,nyngp,nzngp
#ifdef READHDF5  
  Integer(hid_t) :: file_id                               ! File identifier
  Character(len=400) :: filehdf5
  Character(len=16) :: grname
  Character(len=16) :: aname
  Integer(kind=hid_t) ::  gr_id
  integer::nres
#endif

  !------------------------------!
  !READ PARAMETERS AND INITIALIZE!
  !------------------------------!

  !MPI initialization
  call MPI_INIT(mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)


  !!!!!!!ATTENTION IL FAUDRA DISTINGUER ETOILES ET PARTICULES QUAND RUN HYDRO!!!!
  

  !optionally try to improve default params ncut, nxbuffer, nsel if an estimate of an upper bound
  !for npart in halo is known (hardcoded value)
  if(nparthalo_upperbound_estimate_default>0) then
     ncut=ceiling(((nparthalo_upperbound_estimate_default*2.d0/seuil)**(1.d0/3.d0))/2.d0)
     nxbuffer=4*ncut
     nsel=ceiling((2.d0*ncut+1.d0)**3*seuil*4.)
     if(procid==0)then 
        print*,'default params based on an upper bound of npart in halo of',nparthalo_upperbound_estimate_default
        print*,'ncut=',ncut,' nsel=',nsel,' nxbuffer=',nxbuffer
     endif
  endif

  !Read params
  call read_params
  
  !optionally optimize params ncut, nxbuffer, nsel if an estimate of an upper bound
  !for npart in halo is known (user defined value)
  if(nparthalo_upperbound_estimate>0) then
     ncut=ceiling(((nparthalo_upperbound_estimate*2.d0/seuil)**(1.d0/3.d0))/2.d0)
     nxbuffer=4*ncut
     nsel=ceiling((2.d0*ncut+1.d0)**3*seuil*4.)
     if(procid==0)then 
        print*,'used params based on an upper bound of npart in halo of',nparthalo_upperbound_estimate
        print*,'ncut=',ncut,' nsel=',nsel,' nxbuffer=',nxbuffer
     endif
  endif

  if(procID==0)then
     print*,' '
     print*,'Will use following parameters'
     print*,'rep         = ',trim(repository)
     print*,'nx          = ',nx
     print*,'overdensity = ',seuil
     print*,'overdensity2= ',seuil2
     print*,'nmin        = ',nmin
     print*,'ncut        = ',ncut
     print*,'nsel        = ',nsel
     print*,'nxbuffer    = ',nxbuffer, ' (if 0 will use nxzone/2)'
     print*,'nchar       = ',nchar
     print*,' '
  endif
  
  
  !number of dimensions and grid size
  ndim=3
  ny=nx
  nz=nx
  if(procID==0)write(*,*)'Working mesh=',nx,ny,nz
  
#ifndef READHDF5 
     !binary
     !Read headers to get number of cube to read and total number of part

     ! Count the number of cube files present
     ncpu=0
     found=.true.
     do while(found.eqv..true.)
        write(idpchar,'(I5.5)') ncpu
        filenametmp=trim(root)//'_cube_'//idpchar
        inquire(file=filenametmp, exist=found)
        ncpu=ncpu+1   
     enddo

     ncpu=ncpu-1
     if(procID==0)write(*,*)'number of files to read', ncpu
  
  if(ncpu.lt.procNB)then
     print*,'number of mpi tasks should be smaller or equal to number of files to read',ncpu,procNB
     print*,'will stop'
     call MPI_FINALIZE(mpierr)
     stop
  endif

  ! Count the total number of particles in the simulation 
  nthalo=0
  if(procID==0)then
     do kk=0,ncpu-1
        write(idpchar,'(I5.5)') kk
        filenametmp=trim(root)//'_cube_'//idpchar
        open(40,file=filenametmp, status='Old', form='Unformatted')
        read(40) nloc
        close(40)
        nthalo=nloc+nthalo
     enddo
  end if
  call MPI_BCAST( nthalo, 1, MPI_INT, 0,MPI_COMM_WORLD, mpierr )
  nparttotal=nthalo
  if(prociD==0)print*,'nparttotal=',nparttotal

#else
   !HDF5
  filename=trim(root)
  ncpu=procNB  !force ncpu to be equal to procNB

  Call hdf5_init()
  idpchar='00000'
  filehdf5 = trim(root)//'_'//idpchar//'.h5'
  if (procID==0)print*,trim(filehdf5)
  Call hdf5_open_file(filehdf5, file_id)
  
  ! open the root group
  grname = '/'
  Call hdf5_open_group(file_id,grname, gr_id)

  aname='nres'
  Call hdf5_read_attr(gr_id,aname,nres)

  Call hdf5_close_group(gr_id)

  Call hdf5_close_file(file_id)


   nparttotal=int(nres,kind=8)**3 !WARNING TO BE CHANGED IN INTEGER8 for large run


#endif


  !Eeach mpi task will analyze a different set of cubes. 
  !Cubes are analyzed one by one by each mpi task-> could be improved
  !Advice to be more efficient: use a number of task equal to number of cube
  !firstp=first cube that will be treated
  !lastp = last one that will be treated
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

! Loop on the utile zones (cube of number nzone) that will be analyzed 
  do nzone=firstp,lastp
     
     !-------------------------------------------------------------------------------!
     !PARALLEL READING OF CUBE FILES AND SELECTION OF PARTICLES (WORKING+BUFFER ZONE)!
     !-------------------------------------------------------------------------------!
     
     ! Move to a 3D index
     nxcpu=idnint(real(ncpu,kind=8)**(1.d0/3.d0))
     nycpu=idnint(real(ncpu,kind=8)**(1.d0/3.d0))
     nzcpu=idnint(real(ncpu,kind=8)**(1.d0/3.d0)) 
     
     ixzone=int(nzone/(real(nycpu*nzcpu,kind=8)))
     iyzone=int((nzone-ixzone*nycpu*nzcpu)/real(nzcpu,kind=8))
     izzone=nzone-ixzone*nycpu*nzcpu-iyzone*nzcpu
     
     !number of cells in one zone along x direction
     nxzone=nx/nxcpu

     !range of index of cubes to read
     ixmin=ixzone-1
     ixmax=ixzone+1
     iymin=iyzone-1
     iymax=iyzone+1
     izmin=izzone-1
     izmax=izzone+1

     if(procID==0)write(*,*) 'Proc', procID, ' Reading utile cube file: ',nzone,' while other reading others'
     if(procID==0)write(*,*) 'Proc', procID, ' 3D indexes',ixzone,iyzone,izzone

     !Compute cube list to read
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
              list_cpu(kk)=mod(ktmp,nzcpu)+mod(jtmp,nycpu)*nzcpu+mod(itmp,nxcpu)*nzcpu*nycpu
              kk=kk+1
           enddo
        enddo
     enddo

     !if not set the buffer zone is half the length of the working zone
     if(nxbuffer==0)then 
        nxbuffer=nxzone/2  
        if(procID==0)print*,'nxbuffer set to half the buffer zone=',nxbuffer
     endif

     !Check validity of buffer zone
     if(nxbuffer < ncut)then
        print*,'be careful nxbuffer should at least be larger or equal than ncut',nxbuffer,ncut
        print*,'will stop'
        call MPI_FINALIZE(mpierr)
        stop
     endif
     if(nxbuffer>nxzone)then
        print*,'be careful nxbuffer should be smaller than nxzone',nxbuffer,nxzone
        print*,'will stop: TODO: handle this case'
        call MPI_FINALIZE(mpierr)
        stop
     endif

     
     !computation for the limits in xmin,xmax of the entire zone studied + center
     xmin=(ixzone-real(nxbuffer,kind=8)/nxzone)/real(nxcpu,kind=8)
     xmax=(ixzone+1+real(nxbuffer,kind=8)/nxzone)/real(nxcpu,kind=8)
     ymin=(iyzone-real(nxbuffer,kind=8)/nxzone)/real(nycpu,kind=8)
     ymax=(iyzone+1+real(nxbuffer,kind=8)/nxzone)/real(nycpu,kind=8)
     zmin=(izzone-real(nxbuffer,kind=8)/nxzone)/real(nzcpu,kind=8)
     zmax=(izzone+1+real(nxbuffer,kind=8)/nxzone)/real(nzcpu,kind=8)

     xmincenter=(ixzone-real(nxbuffer,kind=8)/(2*nxzone))/real(nxcpu,kind=8)
     xmaxcenter=(ixzone+1+real(nxbuffer,kind=8)/(2*nxzone))/real(nxcpu,kind=8)
     ymincenter=(iyzone-real(nxbuffer,kind=8)/(2*nxzone))/real(nycpu,kind=8)
     ymaxcenter=(iyzone+1+real(nxbuffer,kind=8)/(2*nxzone))/real(nycpu,kind=8)
     zmincenter=(izzone-real(nxbuffer,kind=8)/(2*nxzone))/real(nzcpu,kind=8)
     zmaxcenter=(izzone+1+real(nxbuffer,kind=8)/(2*nxzone))/real(nzcpu,kind=8)

     !Boundary conditions: these quantities can be below 0 or above 1

     if(procID==0)write(*,*) 'Proc', procID, ' Zone (including buffer):', xmin,xmax,ymin,ymax,zmin,zmax

     !limits in xmin,xmax without the buffer zone
     xminzone=ixzone/real(nxcpu,kind=8)
     xmaxzone=(ixzone+1)/real(nxcpu,kind=8)
     yminzone=iyzone/real(nycpu,kind=8)
     ymaxzone=(iyzone+1)/real(nycpu,kind=8)
     zminzone=izzone/real(nzcpu,kind=8)
     zmaxzone=(izzone+1)/real(nzcpu,kind=8)

     !center of the zone
     xc=(xminzone+xmaxzone)/2.d0
     yc=(yminzone+ymaxzone)/2.d0
     zc=(zminzone+zmaxzone)/2.d0

     !distances between the center and the limits
     rxsmall=(0.5d0+real(nxbuffer,kind=8)/real(nxzone,kind=8))/real(nxcpu,kind=8)
     rysmall=(0.5d0+real(nxbuffer,kind=8)/real(nxzone,kind=8))/real(nycpu,kind=8)
     rzsmall=(0.5d0+real(nxbuffer,kind=8)/real(nxzone,kind=8))/real(nzcpu,kind=8)

     rxlarge=1.d0-rxsmall
     rylarge=1.d0-rysmall
     rzlarge=1.d0-rzsmall

     if(procID==0)write(*,*)'Proc', procID, ' Working zone (no buffer):', xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone

     !count number of particles in the selected zone
     
     ! Wait for the token                                                                                                                                                                              
     if(IOGROUPSIZE>0) then
        if (mod(procID,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,procID-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
        end if
     endif


#ifndef READHDF5
        !START READING BINARY
        if (prociD==0)print*,'READ BINARY'
        npart=0
        countinfile: DO kk=0,26
           !read pos
           write(idpchar,'(I5.5)') list_cpu(kk)
           filenametmp=trim(root)//'_cube_'//idpchar        
           open(uinout,file=filenametmp, status='Old', form='Unformatted')
           read(uinout) nloc
           read(uinout) proc_ID
           read(uinout)mintmp(1),maxtmp(1),mintmp(2),maxtmp(2),mintmp(3),maxtmp(3)
           allocate(cdm(3,nloc))
           read(uinout) ((cdm(j,i),j=1,3),i=1,nloc)
           close(uinout)
           
           !count number of particles
           do i=1,nloc
              dx=abs(cdm(1,i)-xc)
              dy=abs(cdm(2,i)-yc)
              dz=abs(cdm(3,i)-zc)
              if ((dx.lt.rxsmall.or.dx.gt.rxlarge).and.&
                   (dy.lt.rysmall.or.dy.gt.rylarge).and.&
                   (dz.lt.rzsmall.or.dz.gt.rzlarge)) then
                 npart=npart+1
              endif
           end do
           deallocate(cdm)
           
        end do countinfile
        
        !allocate final arrays
        allocate(posf(1:ndim,1:npart),velf(1:ndim,1:npart))
        allocate(idf(1:npart))
        allocate(mpf(1:npart))
        
        !fill final arrays x, vel, id, mp
        ifillpart=0
        loadfile: DO kk=0,26
           !read
           write(idpchar,'(I5.5)') list_cpu(kk)
           filenametmp=trim(root)//'_cube_'//idpchar
           open(uinout,file=filenametmp, status='Old', form='Unformatted')
           read(uinout) nloc
           read(uinout) proc_ID
           read(uinout)mintmp(1),maxtmp(1),mintmp(2),maxtmp(2),mintmp(3),maxtmp(3)
           
           allocate(cdm(3,nloc),v(3,nloc),idtmp(nloc))
           read(uinout) ((cdm(j,i),j=1,3),i=1,nloc)
           read(uinout) ((v(j,i),j=1,3),i=1,nloc)
           read(uinout) (idtmp(i),i=1,nloc)
           close(uinout)
           
           
           !fill
           do i=1,nloc
              dx=abs(cdm(1,i)-xc)
              dy=abs(cdm(2,i)-yc)
              dz=abs(cdm(3,i)-zc)
              if ((dx.lt.rxsmall.or.dx.gt.rxlarge).and.&
                   (dy.lt.rysmall.or.dy.gt.rylarge).and.&
                   (dz.lt.rzsmall.or.dz.gt.rzlarge)) then
                 ifillpart=ifillpart+1
                 
                 do j=1,3
                    posf(j,ifillpart)=cdm(j,i)
                    velf(j,ifillpart)=v  (j,i)
                 enddo
                 mpf      (ifillpart)=1.  !non zoom sims
                 idf(ifillpart)=idtmp(i)
                 
              endif
              
           end do
           
           !deallocate temporary arrays
           deallocate(cdm,v,idtmp)
        ENDDO loadfile
#else
      
        !read HDF5 to get npart,pos,vel,mp,id
        if (prociD==0)print*,'READ HDF5'        
        type = 'cuboid'

        center(1)=xc
        center(2)=yc
        center(3)=zc
        dimensions(1)=xmax-xmin
        dimensions(2)=ymax-ymin
        dimensions(3)=zmax-zmin
       
        if(procID==0) then
           print*,type
           print*,filename
           print*,center
           print*,dimensions
        endif
        Call extract_cuboid()
        
        npartini=size(id)

        if(procID ==0)print*,procID,'npart read',npartini
        
        
        !Select only particles in region=> TODO to be improved for memory
        !count number of particles
        npart=0
        do i=1,npartini
           if(procID==0) then
              write(9999,*)pos(1,i),pos(2,i),pos(3,i)
           endif
           dx=abs(pos(1,i)-xc)
           dy=abs(pos(2,i)-yc)
           dz=abs(pos(3,i)-zc)
           if ((dx.lt.rxsmall.or.dx.gt.rxlarge).and.&
                (dy.lt.rysmall.or.dy.gt.rylarge).and.&
                (dz.lt.rzsmall.or.dz.gt.rzlarge)) then
              npart=npart+1
           endif
        end do
        

        if(procID==0)then
           print*,procID, 'npart count region=',npart
           print*,'rxsmall',rxsmall
           print*,'rxlarge',rxlarge
           print*,'center',xc,yc,zc
        endif
        
        allocate(posf(1:ndim,1:npart),velf(1:ndim,1:npart))
        allocate(idf(1:npart))
        allocate(mpf(1:npart))
        
        ifillpart=0
        do i=1,npartini
           dx=abs(pos(1,i)-xc)
           dy=abs(pos(2,i)-yc)
           dz=abs(pos(3,i)-zc)
           if ((dx.lt.rxsmall.or.dx.gt.rxlarge).and.&
                (dy.lt.rysmall.or.dy.gt.rylarge).and.&
                (dz.lt.rzsmall.or.dz.gt.rzlarge)) then
              ifillpart=ifillpart+1
              
              do j=1,3
                 posf(j,ifillpart)=pos(j,i)
                 velf(j,ifillpart)=vel(j,i)
              enddo
              mpf      (ifillpart)=1.  !non zoom sims
              idf(ifillpart)=id(i)
              
           endif           
        end do
        deallocate(pos,vel,id)
        
#endif

        if(procID==0)print*,'Proc',procID,' npart=',npart

     ! Send the token                                                                                                                                                                              
     if(IOGROUPSIZE>0) then
        if(mod(procID+1,IOGROUPSIZE)/=0 .and.((procID+1).lt.procNB))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,procID+1,tag, &
                & MPI_COMM_WORLD,info)
        end if
     endif


     !END OF READING


     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
     
     if(procID==0)then
        write(*,*) 'For all proc: Particles have been read + selection particles working/buffer zone'
        write(*,*) 'Proc', procID, 'Found ', npart,' particles in the working+buffer zone (periodic BC)'
        print*,'Proc',procid,'xmin xmax',xmin,xmax
        print*,'Proc',procid,'minmax(x)',minval(posf(1,:)),maxval(posf(1,:))
        print*,'Proc',procid,'ymin ymax',ymin,ymax
        print*,'Proc',procid,'minmax(y)',minval(posf(2,:)),maxval(posf(2,:))
        print*,'Proc',procid,'zmin zmax',zmin,zmax
        print*,'Proc',procid,'minmax(z)',minval(posf(3,:)),maxval(posf(3,:))
     endif
     
     
     !------------------------------------------------!
     !COMPUTATION OF THE DENSITY OF PARTICLES WITH CIC!
     !------------------------------------------------!
     allocate(xsort(1:npart),isort(1:npart))
     !Call CIC (see parameters meaning inside routine)
     if (docic)then
        mpf=mpf/real(nparttotal,kind=8) !conversion forward could be improved
        
        !cic parameters
        xxmin=xmin;  xxmax=xmax;  yymin=ymin ;yymax=ymax ; zzmin=zmin ; zzmax=zmax
        xmincic=xmin;  xmaxcic=xmax ;ymincic=ymin ; ymaxcic=ymax ; zmincic=zmin; zmaxcic=zmax
        nxcic=2*nx; nycic=2*ny ; nzcic=2*nz
        local_nx=2*(nxzone+2*nxbuffer);  local_ny=2*(nxzone+2*nxbuffer);  local_nz=2*(nxzone+2*nxbuffer)
        ixzonecic=ixzone;iyzonecic=iyzone;izzonecic=izzone
        nxzonecic=nxcpu;nyzonecic=nycpu;nzzonecic=nzcpu
        periodic=.true.

        !call cic to compute density in cube
        call cic(npart,xxmin,xxmax,yymin,yymax,zzmin,zzmax,xmincic,xmaxcic,ymincic,ymaxcic,zmincic,zmaxcic,&
             &local_nx,local_ny,local_nz,nxcic,nycic,nzcic,ixzonecic,iyzonecic,izzonecic,nxzonecic,nyzonecic,nzzonecic,periodic)
        if(procID==0)print*,procid,'minmeanmax(cube)',minval(cube),sum(cube)/(1.d0*(local_nx+1)*(local_ny+1)*(local_nz+1)),maxval(cube)

        !call cic^-1 to compute particles density in rho
        call cicm1(npart,xxmin,xxmax,yymin,yymax,zzmin,zzmax,local_nx,local_ny,local_nz,nxcic,nycic,nzcic,periodic)
        if(procID==0)print*,'minmeanmax(rho)',minval(rho),sum(rho)/(1.d0*npart),maxval(rho)
        
        !deallocate
        deallocate(cube)
        
        mpf=mpf*real(nparttotal,kind=8) !conversion backward could be improved
     endif
     
     !conversion for next routines
     scale=real(nx,kind=8)
     posf=posf*scale
     xminzone=xminzone*scale
     xmaxzone=xmaxzone*scale
     yminzone=yminzone*scale
     ymaxzone=ymaxzone*scale
     zminzone=zminzone*scale
     zmaxzone=zmaxzone*scale
     !this is just a factor 1.
     Mmin=dble(nmin)

     !------------------------------------!
     ! SORT PARTICLES ACCORDING TO DENSITY!
     !------------------------------------!
     if(docic)then
        if(procID==0)write(*,*) 'Sorting particles by CIC density'
        call tri_cic(xsort,isort,npart,nxcic,nycic,nzcic,npart_new,seuil)
     else
        if(procID==0)write(*,*) 'From', procID, 'Sorting particles by NGP density'
        nxngp=4*nx
        nyngp=4*ny
        nzngp=4*nz
        call tri_ngp(xsort,isort,npart,nxngp,nyngp,nzngp,npart_new)
     endif
     
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
     
     if(procID==0)then
        write(*,*) 'For all proc: Particles densities have been computed and particles are sorted'
     endif

     !---------------------------------------------------------
     ! HALO DETECTION USING SPHERICAL OVERDENSITY
     !---------------------------------------------------------
     if(procID==0)then
        write(*,*) 'Detect halo'
        write(*,*) 'using a spherical overdensity delta_bar=',seuil
        write(*,*) 'and measure mass also for a second overdensity delta_bar=',seuil2
     end if

     call spherover(nchar,npart,nx,ny,nz,&
          npart_new,ncut,seuil,seuil2,Mmin,scale,nparttotal,nzone,&
          xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone,&
          xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter,nsel)

     deallocate(posf,velf,idf)
     deallocate(mpf)
     deallocate(xsort,isort)
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

     if(procID==0)then
        write(*,*) 'For all proc: this working zone is analyzed'
        write(*,*) 'Ex: proc',procID,'started',firstp,'now at',nzone,'end is',lastp
     endif
     
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  if(procID==0)write(*,*) 'Run completed!'

  Call hdf5_finalize() 
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

      !read parameter on the line command
      n = iargc()
      if (n < 6) then
         if(prociD==0) then
         print *, 'usage: sod  [-inp input_dir (input directory with cube files)]'
         print *, '            [-nx  nx_grid   (number of grid for selection, cic will be twice this)]  '
         print *, '            [-nch nchar for output files of the form sod_nchar*]  '
         print *, '            [-del delta]  (optional: overdensity threshold)'
         print *, '            [-dl2 delta2] (optional: second threshold to evaluate mass inside)'
         print *, '            [-min np_min] (optional: minimum number of particles in halo)'
         print *, ' '
         print *, 'Then if you want you can specify some of these 3 buffer parameters'
         print *, '            [-sel ncut]     (optional:maximum half extension of halo in number of grid)'
         print *, '            [-nse nsel]     (optional:maximum number of particles within cube of size 2*ncut+1)'
         print *, '            [-nxb nxbuffer] (optional:buffer zone for mpi: several times ncut )'
         print *, ' or an upper bound of the number of particles in the most massive halo'
         print *, '            [-npu ] nparthalo_upperbound_estimate (optional: upper bound for halo nb part)'
         print *, ' if you specify both highest priority is given to the upper bound of npart'
         print *, ' '
         print *, '            [-hlp]          (optional)'
         endif
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
         case('-nse')
            read (arg,*) nsel
         case ('-del') 
            read (arg,*) seuil
         case ('-dl2')
            read (arg,*) seuil2
         case ('-nch')
            read (arg,*) nchar
         case('-nxb')
            read (arg,*) nxbuffer
         case('-npu')
            read (arg,*)nparthalo_upperbound_estimate
            
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
!SORT PARTICLES ACCORDING TO NGP DENSITY (original version)
!------------------------------------------------------------
subroutine tri_ngp(pdens,isort,npart,nx,ny,nz,npart_new)
  use modvar,only:posf,mpf

  use mpi
  implicit none
  integer::mpierr,procID,procNB
  integer::npart,nx,ny,nz,npart_new
  integer,dimension(1:npart)::isort
  real(8),dimension(1:npart)::pdens

  integer(8),dimension(:),allocatable::indx
  integer::ngxngy,ngx,ngy,ngz
  integer(8)::nmax,i1,i2,i3,ind
  integer::i,j,ifail
  integer::imin,imax,ntot
  real(8)::amass
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)
  
  ngx=nx
  ngy=ny
  ngz=nz

  ngxngy=ngx*ngy

  write(*,*) 'npart here is ',npart
   
  allocate(indx(1:npart))

  ! Remplissage du tableau des indices
  !-----------------------------------
  do i = 1,npart
     i1 = int(posf(1,i)*dble(ngx)/dble(nx))
     i2 = int(posf(2,i)*dble(ngy)/dble(ny))
     i3 = int(posf(3,i)*dble(ngz)/dble(nz))
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
     amass= 0.0d0
     do i = imin ,imax
       amass=amass+mpf(isort(i))
     end do
     do i = imin ,imax
        pdens(isort(i)) = 1.0d0/amass
     end do
     imin = imax + 1
  end do

  ! Tri des particules selon la densite
  !------------------------------------
  call quick_sort(pdens,isort,npart)

  i=1
  do while (pdens(i) < 0.5d0)
     i=i+1
  end do
  npart_new=i
  if(procID==0)write(*,*)'npart_active=',npart_new,'/',npart
  deallocate(indx)
  return

end subroutine tri_ngp


subroutine tri_cic(pdens,isort,npart,nx,ny,nz,npart_new,seuil)
  use mpi
  use modvar,only:posf,mpf
  use modcic,only:rho
  implicit none
  integer::npart,nx,ny,nz,npart_new
  integer,dimension(1:npart)::isort
  real(8),dimension(1:npart)::pdens
  real(8)::seuil
  
  integer(8),dimension(:),allocatable::indx
  !integer(8),dimension(1:npart)::indx
  integer::ngxngy,ngx,ngy,ngz
  integer(8)::nmax,i1,i2,i3,ind
  integer::i,j,ifail
  integer::imin,imax,ntot
  real(8)::amass
  integer :: procNB,procID,mpierr
  call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)


  ! Tri des particules selon la densite
  !------------------------------------
  pdens=1.d0/rho
  deallocate(rho)

  call quick_sort(pdens,isort,npart)
  i=1
  !do while (pdens(i) < 1.d0/2.d0) !the 2 here can be increase up to delta or a fraction of delta I think RY
  do while (pdens(i) < 1.d0/(seuil/2.d0)) !RY
     i=i+1
  end do
  npart_new=i
  if(procID==0)write(*,*)'npart_active=',npart_new,'/',npart

  return

end subroutine tri_cic

subroutine cic(npart,xxmin,xxmax,yymin,yymax,zzmin,zzmax,xmin,xmax,ymin,ymax,zmin,zmax,&
     &local_nx,local_ny,local_nz,nx,ny,nz,ixzone,iyzone,izzone,nxzone,nyzone,nzzone,periodic)
!Compute density cube from position array x and mass mp using CIC on a grid of nx_local,
!ny_local,nz_local in utile zone between xxmin, xxmax, yymin,yymax,zzmin,zzmax
!The full region where particles are considered is xmin,xmax,ymin,ymax,zmin,zmax
!The additionnal buffer zone could be used for parallelisation.
!Units: position between [0.,1.], masses so that the sum of the masses of all  particles
!in the volume [0.,1.]^3 is 1 if the density is the average matter universe density.
!nx,ny,nz is the full grid size that covers the full simulation between 0 and 1.
!ixzone, iyzone, izzone is the position of the utile zone in the full volume 
!(ex: this is the zeroth zone along x, the second along y and the second along z)
!nxzone,nyzone,nzzone is the total number of zone in each direction
!Periodic allows period BC to be included 
!This CIC routine is adapted from powergrid.
!Y.Rasera September, 2015
!TO DO: allows the full volume to be different from [0.,1.] to be more flexible
use mpi
use modvar,only:posf,mpf

                                             !input local array x  (positions between [0.,1.[)
                                             !input local array mp (masses, 1./nparttot for non zoom simulation pure DM cosmological simulations)

use modcic,only:cube                         !output local array cube (cic density)

implicit none
integer::npart                               !input local number of particles
real(8)   ::xxmin,xxmax,yymin,yymax,zzmin,zzmax !input local utile zone
real(8)   ::xmin ,xmax ,ymin ,ymax ,zmin ,zmax  !input local zone with buffer
integer::local_nx,local_ny,local_nz          !input local utile grid size
integer::nx,ny,nz                            !input total grid size
integer::ixzone,iyzone,izzone                !position of utile zone in the full simulation
integer::nxzone,nyzone,nzzone                !total number of utile zone in each direction
logical::periodic                            !input periodic boundary conditions


integer::idim=1,jdim=2,kdim=3                !historical: was direction of projection
real(8)::dx,dy,dz,ddx,ddy,ddz,dex,dey,dez
integer::ix,iy,iz,ixp1,iyp1,izp1
integer::i,ntmp
logical::ok_part
real(8)::mtmp
integer :: procNB,procID,mpierr

call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)

allocate(cube(0:local_nx,0:local_ny,0:local_nz))
cube=0.d0

dx=1.d0/real(nx,kind=8)
dy=1.d0/real(ny,kind=8)
dz=1.d0/real(nz,kind=8)
mtmp=0.d0
ntmp=0

!rem: if xxmin==xmin, the density will be correct except near the edge
if(periodic)then
   !loop on all particles
   do i=1,npart
      if (posf(1,i)==1.)write(*,*)'WARNING x=1'
      if (posf(2,i)==1.)write(*,*)'WARNING y=1'
      if (posf(3,i)==1.)write(*,*)'WARNING z=1'
      !selection of particles of interest
      ok_part=(((posf(1,i)>=xmin   .and.posf(1,i)<=xmax   ).or.     &
           &    (posf(1,i)>=xmin+1.d0.and.posf(1,i)<=xmax+1.d0).or.     &
           &    (posf(1,i)>=xmin-1.d0.and.posf(1,i)<=xmax-1.d0))   .and.&
           &   ((posf(2,i)>=ymin   .and.posf(2,i)<=ymax   ).or.     &
           &    (posf(2,i)>=ymin+1.d0.and.posf(2,i)<=ymax+1.d0).or.     &
           &    (posf(2,i)>=ymin-1.d0.and.posf(2,i)<=ymax-1.d0))   .and.&
           &   ((posf(3,i)>=zmin   .and.posf(3,i)<=zmax   ).or.     &
           &    (posf(3,i)>=zmin+1.d0.and.posf(3,i)<=zmax+1.d0).or.     &
           &    (posf(3,i)>=zmin-1.d0.and.posf(3,i)<=zmax-1.d0)))
      
      if(ok_part)then
         ntmp=ntmp+1
         mtmp=mtmp+mpf(i)
         !compute distances from edges of cells in all direction
         ddx=(posf(idim,i)-xxmin)/dx
         ddy=(posf(jdim,i)-yymin)/dy
         ddz=(posf(kdim,i)-zzmin)/dz
         ix=ddx
         iy=ddy
         iz=ddz
         !We subtract because conversion to integer is different 
         !when going to negative numbers. We want a floor.
         if (ddx < 0.d0)ix=ix-1 
         if (ddy < 0.d0)iy=iy-1
         if (ddz < 0.d0)iz=iz-1
         ddx=ddx-ix
         ddy=ddy-iy
         ddz=ddz-iz
         dex=1.0d0-ddx
         dey=1.0d0-ddy
         dez=1.0d0-ddz
         !compute index of neighbouring cells
         ixp1=ix+1
         iyp1=iy+1
         izp1=iz+1

         !periodic BC (assumes CIC + utile zone is within [0,1[+ serial)
         if(ixp1<0)ixp1=ixp1+nx 
         if(ix<0)ix=ix+nx
         if(ixp1>=nx)ixp1=ixp1-nx 
         if(ix>=nx)ix=ix-nx

         if(iyp1<0)iyp1=iyp1+ny 
         if(iy<0)iy=iy+ny
         if(iyp1>=ny)iyp1=iyp1-ny 
         if(iy>=ny)iy=iy-ny

         if(izp1<0)izp1=izp1+nz 
         if(iz<0)iz=iz+nz
         if(izp1>=nz)izp1=izp1-nz 
         if(iz>=nz)iz=iz-nz
         
         
         !Compute density
         if (iz>=0.and.iz<=local_nz) then 
            if (iy>=0.and.iy<=local_ny) then
               if (ix>=0.and.ix<=local_nx)    cube(ix  ,iy  ,iz  )=cube(ix  ,iy  ,iz  )+mpf(i)*dex*dey*dez
               if (ixp1>=0.and.ixp1<=local_nx)cube(ixp1,iy  ,iz  )=cube(ixp1,iy  ,iz  )+mpf(i)*ddx*dey*dez
            endif
            if (iyp1>=0.and.iyp1<=local_ny) then
               if (ix>=0.and.ix<=local_nx)    cube(ix  ,iyp1,iz  )=cube(ix  ,iyp1,iz  )+mpf(i)*dex*ddy*dez
               if (ixp1>=0.and.ixp1<=local_nx)cube(ixp1,iyp1,iz  )=cube(ixp1,iyp1,iz  )+mpf(i)*ddx*ddy*dez
            endif
         endif
         if (izp1>=0.and.izp1<=local_nz) then
            if (iy>=0.and.iy<=local_ny) then
               if (ix>=0.and.ix<=local_nx)    cube(ix  ,iy  ,izp1)=cube(ix  ,iy  ,izp1)+mpf(i)*dex*dey*ddz
               if (ixp1>=0.and.ixp1<=local_nx)cube(ixp1,iy  ,izp1)=cube(ixp1,iy  ,izp1)+mpf(i)*ddx*dey*ddz
            endif
            if (iyp1>=0.and.iyp1<=local_ny) then
               if (ix>=0.and.ix<=local_nx)    cube(ix  ,iyp1,izp1)=cube(ix  ,iyp1,izp1)+mpf(i)*dex*ddy*ddz
               if (ixp1>=0.and.ixp1<=local_nx)cube(ixp1,iyp1,izp1)=cube(ixp1,iyp1,izp1)+mpf(i)*ddx*ddy*ddz
            endif
         endif
      end if
   end do
else
   do i=1,npart
      ok_part=(posf(1,i)>=xmin.and.posf(1,i)<xmax.and. &
           &   posf(2,i)>=ymin.and.posf(2,i)<ymax.and. &
           &   posf(3,i)>=zmin.and.posf(3,i)<zmax)
      if(ok_part)then
         !This section computes distance to edges of cells
         ddx=(posf(idim,i)-xxmin)/dx
         ddy=(posf(jdim,i)-yymin)/dy
         ddz=(posf(kdim,i)-zzmin)/dz
         ix=ddx
         iy=ddy
         iz=ddz
         !because conversion to integer is different when going to negative numbers. 
         !We want a floor.
         if (ddx < 0.d0)ix=ix-1 
         if (ddy < 0.d0)iy=iy-1 
         if (ddz < 0.d0)iz=iz-1 
         ddx=ddx-ix
         ddy=ddy-iy
         ddz=ddz-iz
         dex=1.0d0-ddx
         dey=1.0d0-ddy
         dez=1.0d0-ddz
         
         !compute neighboring cells
         ixp1=ix+1
         iyp1=iy+1
         izp1=iz+1
         
         !compute density
         if(ix>=0.and.ix<=nx.and.iy>=0.and.iy<=ny.and.iz>=0.and.iz<=nz)then 
            cube(ix  ,iy  ,iz  )=cube(ix  ,iy  ,iz  )+mpf(i)*dex*dey*dez
            cube(ix  ,iyp1,iz  )=cube(ix  ,iyp1,iz  )+mpf(i)*dex*ddy*dez
            cube(ixp1,iy  ,iz  )=cube(ixp1,iy  ,iz  )+mpf(i)*ddx*dey*dez
            cube(ixp1,iyp1,iz  )=cube(ixp1,iyp1,iz  )+mpf(i)*ddx*ddy*dez
            cube(ix  ,iy  ,izp1)=cube(ix  ,iy  ,izp1)+mpf(i)*dex*dey*ddz
            cube(ix  ,iyp1,izp1)=cube(ix  ,iyp1,izp1)+mpf(i)*dex*ddy*ddz
            cube(ixp1,iy  ,izp1)=cube(ixp1,iy  ,izp1)+mpf(i)*ddx*dey*ddz
            cube(ixp1,iyp1,izp1)=cube(ixp1,iyp1,izp1)+mpf(i)*ddx*ddy*ddz
         endif
      end if
   end do
endif


!divide by volume
cube=cube/(dx*dy*dz)

end subroutine cic


subroutine cicm1(npart,xxmin,xxmax,yymin,yymax,zzmin,zzmax,local_nx,local_ny,local_nz,nx,ny,nz,periodic)
!Compute particles density rho from position array x and grid density cube using CIC^-1.
!The utile zone is xxmin,xxmax  yymin,yymax zzmin,zzmax 
!It is extracted from a simulation between [0,1.] with grid of size nx,ny,nz
!Units: position between [0.,1.], masses so that the sum of the masses of all  particles
!in the volume [0.,1.]^3 is 1 if the density is the average matter universe density.
!This CIC^-1 routine is adapted from powergrid.
!Y.Rasera September, 2015
use mpi
  use modvar,only:posf                        !input local array x  (positions between [0.,1.[)
use modcic,only:cube,rho                     !input local array cube (cic density)
                                             !output local array rho (cic^-1 density for each part)
implicit none
integer::npart                               !input local number of particles
real(8)   ::xxmin,xxmax,yymin,yymax,zzmin,zzmax !input local utile zone
integer::nx,ny,nz                            !full simulation grid size
integer::local_nx,local_ny,local_nz          !input local utile grid size
logical::periodic                            !input periodic boundary conditions

integer::idim=1,jdim=2,kdim=3                !historical: was direction of projection
real(8)::dx,dy,dz,ddx,ddy,ddz,dex,dey,dez
integer::ix,iy,iz,ixp1,iyp1,izp1
integer::i,ntmp
logical::ok_part
integer :: procNB,procID,mpierr
logical::ok_partsafe
call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)

allocate(rho(npart))
rho=0.d0

dx=1.d0/real(nx,kind=8)
dy=1.d0/real(ny,kind=8)
dz=1.d0/real(nz,kind=8)

!rem: if xxmin==xmin, the density will be correct except near the edge (in a region of size dx)

ntmp=0
do i=1,npart
   if(periodic) then
      ok_part=(((posf(1,i)>=xxmin   .and.posf(1,i)<=xxmax   ).or.     &
           &    (posf(1,i)>=xxmin+1.d0.and.posf(1,i)<=xxmax+1.d0).or.     &
           &    (posf(1,i)>=xxmin-1.d0.and.posf(1,i)<=xxmax-1.d0))   .and.&
           &   ((posf(2,i)>=yymin   .and.posf(2,i)<=yymax   ).or.     &
           &    (posf(2,i)>=yymin+1.d0.and.posf(2,i)<=yymax+1.d0).or.     &
           &    (posf(2,i)>=yymin-1.d0.and.posf(2,i)<=yymax-1.d0))   .and.&
           &   ((posf(3,i)>=zzmin   .and.posf(3,i)<=zzmax   ).or.     &
           &    (posf(3,i)>=zzmin+1.d0.and.posf(3,i)<=zzmax+1.d0).or.     &
           &    (posf(3,i)>=zzmin-1.d0.and.posf(3,i)<=zzmax-1.d0)))
   else
      ok_part=(posf(1,i)>=xxmin.and.posf(1,i)<xxmax.and. &
           &   posf(2,i)>=yymin.and.posf(2,i)<yymax.and. &
           &   posf(3,i)>=zzmin.and.posf(3,i)<zzmax)
   endif
   
   
   if(ok_part)then
      ntmp=ntmp+1
      ddx=(posf(idim,i)-xxmin)/dx
      ddy=(posf(jdim,i)-yymin)/dy
      ddz=(posf(kdim,i)-zzmin)/dz
      ix=ddx
      iy=ddy
      iz=ddz
      if (ddx < 0.d0)ix=ix-1 !because conversion to integer is different when going to negative numbers. We want a floor.                                                                          
      if (ddy < 0.d0)iy=iy-1
      if (ddz < 0.d0)iz=iz-1
      ddx=ddx-ix
      ddy=ddy-iy
      ddz=ddz-iz
      dex=1.0d0-ddx
      dey=1.0d0-ddy
      dez=1.0d0-ddz
      ixp1=ix+1
      iyp1=iy+1
      izp1=iz+1

      
      !periodic BC (assumes CIC + utile zone is within [0,1[+ serial)
      if (periodic) then
         if(ixp1<0)ixp1=ixp1+nx 
         if(ix<0)ix=ix+nx
         if(ixp1>=nx)ixp1=ixp1-nx 
         if(ix>=nx)ix=ix-nx
         
         if(iyp1<0)iyp1=iyp1+ny 
         if(iy<0)iy=iy+ny
         if(iyp1>=ny)iyp1=iyp1-ny 
         if(iy>=ny)iy=iy-ny
         
         if(izp1<0)izp1=izp1+nz 
         if(iz<0)iz=iz+nz
         if(izp1>=nz)izp1=izp1-nz 
         if(iz>=nz)iz=iz-nz
      endif
      
      
      if (iz>=0.and.iz<=local_nz) then
         if (iy>=0.and.iy<=local_ny) then
            if (ix>=0.and.ix<=local_nx)    rho(i)= rho(i)+cube(ix  ,iy  ,iz  )*dex*dey*dez
            if (ixp1>=0.and.ixp1<=local_nx)rho(i)= rho(i)+cube(ixp1,iy  ,iz  )*ddx*dey*dez   
         endif
         if (iyp1>=0.and.iyp1<=local_ny) then
            if (ix>=0.and.ix<=local_nx)     rho(i)= rho(i)+cube(ix  ,iyp1,iz  )*dex*ddy*dez
            if (ixp1>=0.and.ixp1<=local_nx) rho(i)= rho(i)+cube(ixp1,iyp1,iz  )*ddx*ddy*dez
         endif
      endif
      if (izp1>=0.and.izp1<=local_nz) then
         if (iy>=0.and.iy<=local_ny) then
            if (ix>=0.and.ix<=local_nx)    rho(i)= rho(i)+cube(ix  ,iy  ,izp1)*dex*dey*ddz
            if (ixp1>=0.and.ixp1<=local_nx)rho(i)= rho(i)+cube(ixp1,iy  ,izp1)*ddx*dey*ddz
         endif
         if (iyp1>=0.and.iyp1<=local_ny) then
            if (ix>=0.and.ix<=local_nx)     rho(i)= rho(i)+cube(ix  ,iyp1,izp1)*dex*ddy*ddz
            if (ixp1>=0.and.ixp1<=local_nx) rho(i)= rho(i)+cube(ixp1,iyp1,izp1)*ddx*ddy*ddz
         endif
      endif
     
   end if
end do

end subroutine cicm1




!---------------------------------------------------------------------------
subroutine spherover(nchar,npart,nx,ny,nz,&
                     npart_new,ncut,seuil,seuil2,Mmin,scale,nparttotal,nzone,&
                     xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone,&
                     xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter,nsel)
!---------------------------------------------------------------------------
  use mpi
  use modvar
#ifdef READHDF5
  use modvariables
  use modhdf5
  use modwritehalo
#endif  
  implicit none
  integer::mpierr,procID,procNB
  integer,parameter:: nbins=60
  character*5::nchar
  integer::npart,nx,ny,nz,npart_new,nparttotal
  real(8)::seuil,Mmin,scale,seuil2
  integer::ncut,ncut2          
  real(8)::xminzone,xmaxzone,yminzone,ymaxzone,zminzone,zmaxzone  
  real(8)::xmincenter,xmaxcenter,ymincenter,ymaxcenter,zmincenter,zmaxcenter  
  integer::Nsel
  integer::np,nzone
  real(8)::twopi=6.283185307179586d0

  integer,dimension(:),allocatable::structure,allpart
  integer,dimension(:),allocatable::masses
  integer,dimension(:),allocatable::indx_ngp,isort_ngp
  integer(idprec),dimension(:),allocatable::idp
  integer,dimension(:),allocatable::first_ngp,num_ngp
  

  integer,dimension(:),allocatable::indx_sel,isort_sel,irank_sel
  real,dimension(:),allocatable::m_sel
  real,dimension(:,:),allocatable::x_sel,x_sel2,v_sel
  real(8),dimension(:),allocatable::distance
  real(8)::rayon_min,distance_new

  integer::i,j,k,l,npart_sel,ipart,npart_sel2,n_max,n_test
  integer::is,namas,gap,masse,ip,ifail,masse2
  integer(idprec)::iamas
  integer::ic,jc,kc,ii,jj,kk,ind,indc
  integer::i1,i2,i3,imin,imax,ntot
  integer::itest,niter=1
  real(8)::xc,yc,zc,xx,yy,zz,xc0,yc0,zc0,dc0,dmax,dmaxtmp
  real(8)::xcnew,ycnew,zcnew,xxnew,yynew,zznew,dxnew,dynew,dznew
  real(8)::xb,yb,zb
  real(8)::taillex,tailley,taillez,amass
  real(8)::dx,dy,dz,d,d2,rayon,rayon2,volume,overdens,overdens2
  real(8)::dmin,x_tmp,y_tmp,z_tmp,dis_tmp,qqq
  real(8)::Mtotsel

  character*5::nzonechar
  character*200::nomf,nomf2,filestrct,fileamas,fileprofile,filestrctfof,filemasst,nomfst
  integer::count_no_halo
  real(8)::frac=1. ! Control when to stop the loop on particles: 1 is safer

  real(8),dimension(1:nbins)::rad_bins
  integer,dimension(1:nbins)::mass_bins
  logical::search_max_npart_in_sphere=.true.
  real(4),dimension(:),allocatable::xctab,yctab,zctab
  integer(idprec),dimension(:),allocatable::idtab
  integer::nmaxhalo
  logical::write_sodtxt=.false.   
  logical::write_profile=.false.  !TO DO: make it working with io ticket
  logical::write_structure=.false.
  logical::write_masst_strct=.true.
  Character(len=400) :: outroot !< root of output 
  
  Integer(idprec), dimension(:), allocatable :: haloID !< Array containing the ID of the halos
  Integer(idprec), dimension(:), allocatable :: halopartID !< Array containing the ID of the
  !< particles belonging to the halos
  Integer(kind=4), dimension(:),allocatable :: haloMass !< Array containing the mass of the halos
  Real   (kind=4), dimension(:,:), allocatable :: halopartPos, halopartVel !< Array containing the position
    !< and velocity of the particles belonging to the halos
  Real(kind=8), dimension(:,:),allocatable:: halocomPos !< Position of the center of mass of the halos
  Real(kind=8), dimension(:,:),allocatable :: halocomVel !< Velocity of the center of mass of the halos
  Real(kind=8), dimension(:),allocatable:: haloRadius !< Radius of the halos
  Integer(kind=4), dimension(:),allocatable:: haloSubHaloNB !< Number of subhalos in each halo
  integer::haloNB_all,nh
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)

  !allocate npart arrays
  allocate(structure(1:npart),indx_ngp(1:npart),isort_ngp(1:npart))
  allocate(allpart(1:npart))
  ! allpart contains the index of the parts of each halo (in x, vel and id arrays) 
  !(from 1 to masses(1), ids of part in first amas, etc)
  
  
  !allocate array with halo properties (for masst files)
  ! masses contains the masses of all the halos found. 
  nmaxhalo=int(dble(npart)/mmin)
  allocate(idtab(1:nmaxhalo))
  allocate(masses(1:nmaxhalo))
  allocate(xctab(1:nmaxhalo))
  allocate(yctab(1:nmaxhalo))
  allocate(zctab(1:nmaxhalo))


  !initialisations
  structure = 0
  taillex   = dble(nx)
  tailley   = dble(ny)
  taillez   = dble(nz)
  rayon_min = sqrt(3.d0)*scale/(2.d0*nx)

  !filenames
  write(nzonechar,'(I5.5)') nzone
  filestrctfof = 'sod_'//trim(nchar)//'_strct_'//trim(nzonechar)
  filemasst    = 'sod_'//trim(nchar)//'_masst_'//trim(nzonechar)

  nomf      = 'sod_'//TRIM(nchar)//'_'//TRIM(nzonechar)//'.dat'
  nomf2     = 'mass_sod_'//TRIM(nchar)//'_'//TRIM(nzonechar)//'.dat'
  fileprofile = 'halos_rad_part_'//trim(nchar)//'_'//TRIM(nzonechar)//'.dat'
  nomfst = 'st_sod_'//TRIM(nchar)//'_'//TRIM(nzonechar)//'.dat'
  
  outroot='sod_'//TRIM(nchar)
  
  !Allocate grid or local part arrays
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

  ! Building particle linked list
  !--------------------------------------------
  if(procID==0)write(*,*)'Building particle linked list...'
  do i  = 1,npart
     i1 = int(posf(1,i))
     i2 = int(posf(2,i))
     i3 = int(posf(3,i))
     indx_ngp(i) = 1+i1+i2*nx+i3*nx*ny
  end do

  xsort=indx_ngp
  call quick_sort(xsort,isort_ngp,npart)
  indx_ngp=xsort

  ! As far as I understand
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

  !print on screen
  if(procID==0)write(*,*)'proc', procID,' found'
  if(procID==0)write(*,*) '  halo#    npart   npart2   mass       x          y          z         it  radius     del   del2'
  
  !open files
  if(write_sodtxt)open(10,file=trim(nomf),form='formatted')
  if(write_sodtxt)open(20,file=trim(nomf2),form='formatted')
  if(write_profile)open(25,file=trim(fileprofile),form='formatted')


  !initialization
  namas = 0
  gap=0
  xsort = 1.0d0
  count_no_halo=0

  ipart=1
  !Loop on particles
  do while( (ipart.le.npart_new).and.(dble(count_no_halo)/dble(npart_new).le.frac))
     
     ! We choose most dense particles
     !-------------------------------------------------
     if(xsort(isort(ipart)).gt.0.)then !this is to eliminate particles which are already in a halo
        
        !We adopt the center of the halo as the current particles
        xc = posf(1,isort(ipart))
        yc = posf(2,isort(ipart))
        zc = posf(3,isort(ipart))

        ! the halo identifier is the id of the barycentre
        iamas = idf(isort(ipart))


        xc0=xc; yc0=yc; zc0=zc !backup the center in case of iteration on the center (old versions)
        
        ! Preselect particles near the center
        !-----------------------------------------------------------
        is = 0
        ic = int(xc)
        jc = int(yc)
        kc = int(zc)
        indc=1+ic+jc*nx+kc*nx*ny
        mtotsel=0.d0
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
                          write(*,*)procID,'WARNING: from ',nzonechar,' Increase Nsel in main program'
                          deallocate(structure,indx_ngp,isort_ngp)
                          deallocate(masses,allpart,idtab,xctab,yctab,zctab)
                          deallocate(num_ngp,first_ngp,indx_sel)
                          deallocate(isort_sel,distance,x_sel,v_sel)
                          deallocate(x_sel2,m_sel,idp)
                          return
                       end if
                       indx_sel(is) =  isort_ngp(ip)
                       x_sel(is,1) = posf(1,isort_ngp(ip))
                       x_sel(is,2) = posf(2,isort_ngp(ip))
                       x_sel(is,3) = posf(3,isort_ngp(ip))
                       v_sel(is,1) = velf(1,isort_ngp(ip))
                       v_sel(is,2) = velf(2,isort_ngp(ip))
                       v_sel(is,3) = velf(3,isort_ngp(ip))
                       m_sel(is) = mpf (isort_ngp(ip))
                       idp(is) = idf(isort_ngp(ip)) !index of the particule
                       mtotsel = mtotsel + m_sel(is)
                    end if
                 end do
              end do
           end do
        end do
        npart_sel = is

        !select particules that are in this cell plus 1 cell around for the barycentre loop 
        !as prospective center
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
                       x_sel2(is,1) = posf(1,isort_ngp(ip))
                       x_sel2(is,2) = posf(2,isort_ngp(ip))
                       x_sel2(is,3) = posf(3,isort_ngp(ip))
                    end if
                 end do
              end do
           end do
        end do
        npart_sel2 = is



        if(Mtotsel .ge. Mmin) then

           !loop on neighbouring particles to find the one with more neighbour
           !computed in the same way as one computes profile=> center is max density

           !compute distance^2 dc0
           dx = abs(xc-xc0)
           dx = min(dx,taillex-dx)
           dy = abs(yc-yc0)
           dy = min(dy,tailley-dy)
           dz = abs(zc-zc0)
           dz = min(dz,taillez-dz)
           dc0 = dx*dx+dy*dy+dz*dz

           if(search_max_npart_in_sphere) then
              
              !---------------------------------------------------
              !loop on all particules in the same cell as the current barycentre
              dmax=1d30
              n_max=0
              l=first_ngp(indc)
              do while(indx_ngp(l).eq.indc.and.l.lt.npart)
                 xcnew=posf(1,isort_ngp(l))
                 ycnew=posf(2,isort_ngp(l))
                 zcnew=posf(3,isort_ngp(l))
                 
                 !compute the distances to the new barycentre for the preselected particules
                 ! and count the number of particules in the sphere of radius rayon_min
                 i=0
                 dmaxtmp=0.d0
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
                    if (distance_new.le.rayon_min**2) then
                       i=i+1
                       if(distance_new.ge.dmaxtmp**2)dmaxtmp=sqrt(distance_new)
                    end if
                 end do
                 n_test=i
                 !update the barycentre if there are more particules in the sphere
                 if (n_test.ge.n_max) then
                    if(n_test.eq.n_max)then   !to get unique answer whatever ncpu
                       if (dmaxtmp.le.dmax) then
                          n_max=n_test
                          xc=xcnew
                          yc=ycnew
                          zc=zcnew
                          iamas=idf(isort_ngp(l))
                          dmax=dmaxtmp
                       endif
                    else
                       n_max=n_test
                       xc=xcnew
                       yc=ycnew
                       zc=zcnew
                       iamas=idf(isort_ngp(l))
                       dmax=dmaxtmp
                    endif
                 end if
                 l=l+1
              end do
           endif
           
           ! Compute distances to the center
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
          
           ! Sort distances
           !---------------------------------------------------
           call quick_sort(distance,isort_sel,npart_sel)

           ! Computation of mass and virial radius
           !---------------------------------------------------
           i        = 1
           overdens = 2.0d0*seuil
           amass    = m_sel(isort_sel(i))
           rayon    = 0.0d0

           do while(i.lt.10.or.(overdens.gt.seuil.and.i.lt.npart_sel.and.rayon.lt.(real(ncut,kind=8)-sqrt(dc0))))
              i        = i+1
              rayon    = SQRT(distance(i))
              volume   = 2.d0/3.d0*twopi*rayon**3
              amass    = amass+m_sel(isort_sel(i))
              overdens = amass/volume
              ! Save mass above a second threshold (seuil2>seuil)
              !---------------------------------------------------------------------------------------
              if (overdens.gt.seuil2) then 
                 masse2=i
                 if (rayon>1d-30) then
                    overdens2=overdens
                 else
                    overdens2=0.d0
                 endif
              end if
           end do

           if(overdens.gt.seuil)then
              write(*,*)procID,'WARNING: from ',nzonechar,' Increase selection radius (-sel option in command line)'
              write(*,*)procID,'mass     =',i
              write(*,*)procID,'nsel     =',npart_sel
              write(*,*)procID,'overdens =',overdens
              write(*,*)procID,'radius   =',rayon+dc0
              write(*,*)procID,'nsel     =',ncut
              write(*,*)procID,'niter    =',niter
              write(*,*)procID,'mass2    =',masse2
              write(*,*)procID,'overdens2=',overdens2
              deallocate(structure,indx_ngp,isort_ngp)
              deallocate(masses,allpart,idtab,xctab,yctab,zctab)
              deallocate(num_ngp,first_ngp,indx_sel,isort_sel)
              deallocate(distance,x_sel,v_sel,x_sel2,m_sel,idp)
              return
           else
              rayon =SQRT(distance(i-1))
              volume=2.d0/3.d0*twopi*rayon**3
              amass =amass-m_sel(isort_sel(i))
              if(volume > 0.d0)then
                 overdens=amass/volume
              else
                 overdens=0.d0
              end if
              masse=i-1
           end if

           
           !Additional loop to compute the profile
           do j=1,nbins
              !bins de rayon repartis logarithmiquement entre le rayon et le rayon*0.01
              rad_bins(j)=rayon*exp(log(0.01d0)+(j-1)*(log(1.d0)-log(0.01d0))/(nbins-1.d0))
              if (j.eq.1) then
                 k=1 !because the center is a particle
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
              !save id,npart,x,y,z for masst files
              idtab(namas)=iamas
              masses(namas)=amass
              xctab(namas)=xc/scale
              yctab(namas)=yc/scale
              zctab(namas)=zc/scale
              
              gap=gap+masse 

              !write files
              if(procID==0)write(*,997)iamas,masse,masse2,dble(amass)/scale**3,xc/scale,yc/scale,zc/scale &
                   & ,niter,rayon/scale,overdens,overdens2
              if(write_sodtxt)write(20,998)iamas,masse,masse2,overdens,overdens2

              if(write_profile)write(25,995)iamas,masse,xc/scale,yc/scale,zc/scale
              do i=1,nbins
                 if(write_profile)write(25,994)rad_bins(i)/scale,mass_bins(i)
              end do
              count_no_halo=0
        end if

        ! We remove halo particles from the list 
        !------------------------------------------------------
        ! to do: possible improvement, we do not remove part if too close to the edge of the buffer zone
        ! 
        !------------------------------------------------------
        do i = 1,masse
           if(amass.ge.Mmin.and.xc.le.xmaxzone.and.xc.gt.xminzone&
                           .and.yc.le.ymaxzone.and.yc.gt.yminzone&
                           .and.zc.le.zmaxzone.and.zc.gt.zminzone)then
              structure(indx_sel(isort_sel(i))) = iamas
              allpart(i+gap-masse)= indx_sel(isort_sel(i))
              xsort(indx_sel(isort_sel(i))) = 0.d0
           else !if (xc.le.xmaxcenter.and.xc.gt.xmincenter&
               !.and.yc.le.ymaxcenter.and.yc.gt.ymincenter&
               !.and.zc.le.zmaxcenter.and.zc.gt.zmincenter)then
              structure(indx_sel(isort_sel(i)))=0
              xsort(indx_sel(isort_sel(i))) = 0.d0
           endif
        enddo
	      
        ! Write main halo caracteristics
        !----------------------------------------------------
        if(amass.ge.Mmin.and.xc.le.xmaxzone.and.xc.gt.xminzone.and.&
             yc.le.ymaxzone.and.yc.gt.yminzone.and.&
             zc.le.zmaxzone.and.zc.gt.zminzone)then
           if(write_sodtxt)write(10,999)iamas,masse,amass/scale**3,xc/scale,yc/scale,zc/scale,rayon/scale
        end if
        
     end if
     
  end if
  count_no_halo=count_no_halo+1

  !write progression on screen
  if(npart_new>10)then
     if(mod(ipart,npart_new/10)==0)then   
        if(procID==0)write(*,'(0PF5.1,"% complete")')100.d0*dble(ipart)/dble(npart_new)  
     endif
  endif
  
  ipart=ipart+1
end do
   
  if(write_sodtxt)close(10)
  if(write_sodtxt)close(20)
  if(write_profile)close(25)
  !if no halo for a while can optionnaly stop
  if(dble(count_no_halo)/dble(npart_new).gt.frac)then
     write(*,*)procID,'Seems that there are no more halos... It s however not guaranteed! '
     write(*,*)procID,'Increase the stopping criterion frac=',frac,' if you think there are still some more '
  endif  
  
  
  ! Wait for the token                                                                                                                                                                              
  if(IOGROUPSIZE>0) then
     if (mod(procID,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,procID-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info)
     end if
  endif

  !write structure (old format)
  if(write_structure)then
     open(30,file=trim(nomfst),form='unformatted')
     write(30)structure
     close(30)
  end if
  
  !write masst file (like pFoF)
  if(write_masst_strct)then
     open(27,file=trim(filemasst),form='unformatted')
     write(27) namas
     do i=1,namas      
        write(27)int(idtab(i),kind=8),masses(i),xctab(i),yctab(i),zctab(i)
     end do
     close(27)


     !write strct file (like pFoF)
     posf=posf/scale !conversion forward to be improved
     gap=0
     open(26,file=trim(filestrctfof),form='unformatted')
     write(26) namas
     do i=1,namas 
        gap=gap+masses(i)
        write(26) masses(i)
        write(26) ((posf(k,allpart(j+gap-masses(i))),k=1,3),j=1,masses(i))
        write(26) ((velf(k,allpart(j+gap-masses(i))),k=1,3),j=1,masses(i))
        write(26) (idf(allpart(j+gap-masses(i))),j=1,masses(i))
     enddo
     close(26)
     
     !put the READHDF5 FLAG LATER 
     !TODO: Improve memory usage
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
     
     if(namas>0) then
        allocate(haloMass(namas))
        allocate(haloID(namas))
        allocate(halocomPos(3,namas))
        
        allocate(halocomVel(3,namas))
        halocomVel=0.
        allocate(haloRadius(namas))
        haloRadius=0.
        allocate(haloSubHaloNB(namas))
        haloSubHaloNB=0.
     else
        !Trick if zero halo on some proc
        allocate(haloMass(1))
        allocate(haloID(1))
        allocate(halocomPos(3,1))
        
        allocate(halocomVel(3,1))
        halocomVel=0.
        allocate(haloRadius(1))
        haloRadius=0.
        allocate(haloSubHaloNB(1))
        haloSubHaloNB=0.
     endif
     gap=0
     do i=1,namas 
        gap=gap+masses(i)
        haloMass(i)=masses(i)
        haloID(i)=int(idtab(i),kind=8) 
        halocomPos(1,i)=xctab(i)
        halocomPos(2,i)=yctab(i)
        halocomPos(3,i)=zctab(i)
     enddo

     allocate(halopartPos(3,gap))
     allocate(halopartVel(3,gap))
     allocate(halopartId (  gap))
     
     gap=0
     do i=1,namas 
        gap=gap+masses(i)
        do j=1,masses(i)
           halopartPos(:,j+gap-masses(i))=posf(:,allpart(j+gap-masses(i)))
           halopartVel(:,j+gap-masses(i))=velf(:,allpart(j+gap-masses(i)))
           halopartID (  j+gap-masses(i))=idf (  allpart(j+gap-masses(i)))
        end do
     end do
    if(procID==16)then
        print*,'outroot',trim(outroot)
        print*,'namas',namas
        print*,'gap',gap
        print*,'haloMass',minval(haloMass),maxval(haloMass)
        print*,'haloID  ',minval(haloID),maxval(haloID)
        print*,'haloppos',minval(halopartPos),maxval(halopartPos)
        print*,'halopvel',minval(halopartVel),maxval(halopartVel)
        print*,'halopID ',minval(halopartid) ,maxval(halopartid)
     endif
     

     call h5writehalopart(outroot,MPI_COMM_WORLD, namas, gap, haloMass, haloID, halopartPos, halopartVel, halopartID)



     Call Mpi_AllReduce(namas,haloNB_all,1,MPI_INTEGER,Mpi_Sum,Mpi_Comm_World,mpierr)
      if(procID==0)then 
         print*,'there are ',haloNB_all,' halos'
      endif
      
     
     nh = ubound(haloMass,1)

     if(procID==16)then
        print*,'outroot',trim(outroot)
        print*,'halonbal',haloNB_all
        print*,'namas',namas
        print*,'nh',nh
        print*,'haloMass',minval(haloMass),maxval(haloMass)
        print*,'halocomp',minval(halocompos),maxval(halocompos)
        print*,'halocomv',minval(halocomvel),maxval(halocomvel)
        print*,'halopID ',minval(haloid) ,maxval(haloid)
        print*,'halopR  ',minval(haloRadius) ,maxval(haloRadius)
        print*,'halopSub',minval(halosubhalonb) ,maxval(halosubhalonb)
        
     endif
     
      
     call  mpih5writehalomass(outroot,MPI_COMM_WORLD, haloNB_all, namas,nh, haloMass, halocomPos, halocomVel, &
       haloID, haloRadius, haloSubHaloNB)
    

     deallocate(haloMass)
     deallocate(haloID)
     deallocate(halocomVel)
     deallocate(haloRadius)
     deallocate(haloSubHaloNB)

     deallocate(halopartPos)
     deallocate(halopartVel)
     deallocate(halopartId)
     
     
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
     if(procID==0)print*,'write hdf5 Done'
     posf=posf*scale !conversion backward to be improved
  endif
  
  ! Send the token                                                                                                                                                                              
  if(IOGROUPSIZE>0) then
     if(mod(procID+1,IOGROUPSIZE)/=0 .and.((procID+1).lt.procNB))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,procID+1,tag, &
             & MPI_COMM_WORLD,info)
     end if
  endif
  
  !deallocate
  deallocate(structure,indx_ngp,isort_ngp,num_ngp,first_ngp)
  deallocate(indx_sel,isort_sel,distance,x_sel,x_sel2,m_sel,idp)
  deallocate(v_sel,allpart,masses,idtab,xctab,yctab,zctab)

997 format (I12,1x,I8,1x,I8,4(1x,1pe10.3),1x,I2,1x,6(E10.3,1x))
998 format (I12,2(1x,I8),2(1x,E10.3))
999 format (I12,1x,I8,5(1x,1pe10.3))
996 format (I12,1x,I12)
995 format (I12,1x,I8,3(1x,1pe10.3))
994 format (E10.3,1x,I8)
	
  return
end subroutine spherover

!Routine to convert in format I5.5
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

!Routine to sort with the quicksort method
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
  
  !Routine called by quick sort
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
