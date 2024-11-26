!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================
Module modvariables

  Use modhdf5

#ifdef WITHMPI
  Use mpi
#endif
  
  ! If you do not know the type of integer used for the ID, use h5dump on the first halo file.
  ! h5dump -d haloID test_halo_00000.h5
  ! This will print the dataset haloID containing the ID of each halo in the file. 
  ! You should see something like
  ! HDF5 "test_halo_00000.h5" {
  ! DATASET "haloID" {
  !    DATATYPE  H5T_STD_I64LE
  ! The datatype should contain the number of bits used, i.e. 32 for kind=4 or 64 for kind=8.
#ifdef LONGINT
  Integer, parameter :: PRI = 8
#ifdef WITHMPI
  Integer, parameter :: MPI_PRI = Mpi_Integer8
#endif
#else
  Integer, parameter :: PRI = 4
#ifdef WITHMPI
  Integer, parameter :: MPI_PRI = Mpi_Integer
#endif
#endif

  ! Name of the codes to write as metadata in HDF5 files
  Character(len=16), parameter :: NAME_CONECREATOR_PART='conecreator_part'
  Character(len=16), parameter :: NAME_CONECREATOR_GRAV='conecreator_grav'
#ifdef SOD
  Character(len=16), parameter :: NAME_PFOF_SNAP='psod_snap'
  Character(len=16), parameter :: NAME_PFOF_CONE='psod_cone'
#else
  Character(len=16), parameter :: NAME_PFOF_SNAP='pfof_snap'
  Character(len=16), parameter :: NAME_PFOF_CONE='pfof_cone'
#endif

  !! Variables to handle the HDF5 file
  Integer(kind=hid_t) :: file_id         ! hdf5 id of the opened file
  Integer(kind=hid_t) :: gr_id           ! hdf5 id of the group
  Integer(kind=hid_t) :: prop_id           ! hdf5 id of the properties group
  Integer(kind=hid_t) :: prof_id           ! hdf5 id of the profile group
  
  Integer(kind=hid_t) :: gr_ramses_id
  Integer(kind=hid_t) :: gr_halo_id
  Integer(kind=hid_t) :: gr_pfof_id
  Integer(kind=hid_t) :: gr_input_id
  Integer(kind=hid_t) :: gr_fof_id
  Integer(kind=hid_t) :: gr_output_id

  Character(len=H5STRLEN) :: name              ! name of the dataset/attribute
  Character(len=H5STRLEN) :: groupname         ! name of the group
  Character(len=9)        :: filenbchar        ! string used for the id of the file
  Integer(kind=4)         :: haloNB            ! number of halos in a file
  Integer(kind=4)         :: halopartNB        ! number of particles in a halo
  Integer(kind=PRI), dimension(:), allocatable :: haloID  ! array containing the ID of each halo in the file

  Integer(kind=PRI) :: currenthaloID ! haloID of the halo analyzed in your function
  Integer(kind=4) :: currenthalo     ! index of the halo which varies from 1 to N where N is the total number of haloes you are analyzing
  Integer(kind=4) :: nbhaloanalyzed  ! total number of haloes that you want to analyze 

  Integer(kind=PRI), dimension(:), allocatable :: id        ! id of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: pos          ! position of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: vel          ! velocities of the particles in the halo just read
  Real(kind=4), dimension(:), allocatable :: pot            ! potential of the particles in the halo just read 

  ! Type declaration for information read from cone info file
  Type :: Type_info_cone !< Type for information read from cone info file
     Integer(kind=4) :: ncpu
     Integer(kind=4) :: nstride
     Integer(kind=4) :: nstep_coarse
     Integer(kind=4) :: cone_id
     Integer(kind=4) :: nglobalfile
     Integer(kind=4) :: isfullsky
     Integer(kind=4) :: future
     Real(kind=8) :: aexp
     Real(kind=8) :: observer_x
     Real(kind=8) :: observer_y
     Real(kind=8) :: observer_z
     Real(kind=8) :: observer_rds
     Real(kind=8) :: cone_zlim
     Real(kind=8) :: amax
     Real(kind=8) :: amin
     Real(kind=8) :: zmax
     Real(kind=8) :: zmin
     Real(kind=8) :: dmax
     Real(kind=8) :: dmin
     Real(kind=8) :: dtol
     Real(kind=8) :: thetay
     Real(kind=8) :: thetaz
     Real(kind=8) :: theta
     Real(kind=8) :: phi
     Real(kind=8) :: aendconem2
     Real(kind=8) :: aendconem1
     Real(kind=8) :: aendcone
     Real(kind=8) :: zendconem2
     Real(kind=8) :: zendconem1
     Real(kind=8) :: zendcone
     Real(kind=8) :: dendconem2
     Real(kind=8) :: dendconem1
     Real(kind=8) :: dendcone
  End Type Type_info_cone

  Type, extends(Type_info_cone) :: Type_info_cone_part !< Extension of Type_info_cone for particles cone
     Integer(kind=8) :: npart          !! nglobalcell in the file
     Real(kind=8) :: aexpold
     Real(kind=8) :: zexpold
     Real(kind=8) :: zexp
     Real(kind=8) :: dexpold
     Real(kind=8) :: dexp
  End Type Type_info_cone_part

  Type, extends(Type_info_cone) :: Type_info_cone_grav !< Extension of Type_info_cone for cells cone
     Integer(kind=4) :: nlevel
     Integer(kind=4) :: level_min
     Integer(kind=4) :: level_max
     Integer(kind=8) :: nglobalcell
  End type Type_info_cone_grav

  !! Variables read from the Ramses info file
  Type :: Type_info_ramses           ! Type for information read from Ramses info file
     Integer(kind=4) :: ncpu
     Integer(kind=4) :: ndim         ! number of dimensions
     Integer(kind=4) :: levelmin     ! minimum mesh refinement level
     Integer(kind=4) :: levelmax
     Integer(kind=4) :: ngridmax
     Integer(kind=4) :: nstep_coarse
     Real(kind=8) :: boxlen
     Real(kind=8) :: time
     Real(kind=8) :: aexp
     Real(kind=8) :: h0
     Real(kind=8) :: omega_m
     Real(kind=8) :: omega_l
     Real(kind=8) :: omega_k
     Real(kind=8) :: omega_b
     Real(kind=8) :: unit_l
     Real(kind=8) :: unit_d
     Real(kind=8) :: unit_t
     Real(kind=8) :: delta_vir
  End Type Type_info_ramses

  !! Variables read from pfof input parameters common to snapshot and cone versions
  Type :: Type_parameter_pfof ! Type for pfof input parameters common to snapshot and cone versions
     Character(len=H5STRLEN) :: input_path              ! path to the directory containing the RAMSES files
     Character(len=H5STRLEN) :: part_input_file         ! name of the files containing the particles data
     Character(len=H5STRLEN) :: simulation_name         ! Simulation name to be written in output files
     Integer(kind=4)    :: mmin                         ! minimum mass of the haloes
     Integer(kind=4)    :: mmax                         ! maximum mass of the haloes
     Logical(kind=4)    :: do_read_potential            ! do the RAMSES files contain potential?
     Logical(kind=4)    :: do_read_gravitational_field  ! do the RAMSES files contain force?
     Logical(kind=4)    :: do_unbinding                 ! Unbinding process? (not implemented yet)
     Logical(kind=4)    :: do_subHalo                   ! SubHalo detection? (not implemented yet)
     Logical(kind=4)    :: do_timings                   ! should pFOF perform timings (imply extra synchronizations)
     Real(kind=4)       :: percolation_length           ! FOF percolation length
  End Type Type_parameter_pfof

  Type, extends(Type_parameter_pfof)  :: Type_parameter_pfof_snap ! Extension of Type_parameter_pfof for the snapshot version
     Integer(kind=4)    :: grpsize              ! size of the I/O group used for the RAMSES simulation
     Integer(kind=4)    :: gatherread_factor    ! gather parameter for parallel hdf5 input
     Integer(kind=4)    :: snapshot             ! Ramses output number to be written in output files
     Integer(kind=4)    :: gatherwrite_factor   ! gather parameter for parallel hdf5 output
     Logical(kind=4)    :: do_skip_star         ! do the RAMSES files contain stars?
     Logical(kind=4)    :: do_skip_metal        ! do the RAMSES files contain metalicities?
     Logical(kind=4)    :: do_read_from_cube    ! should pFOF read particles from cube files?
     Logical(kind=4)    :: do_fof               ! should pFOF perform FOF halo finding?
     Logical(kind=4)    :: do_write_cube        ! should pFOF write cube files?
     Logical(kind=4)    :: do_sort_cube         ! sort the particles and write a 'sorted cube'
     Character(len=3)   :: code_index           ! version of RAMSES used to produce the RAMSES files
     Character(len=H5STRLEN) :: info_input_file ! name of the RAMSES info file
  End Type Type_parameter_pfof_snap

  Type, extends(Type_parameter_pfof) :: Type_parameter_pfof_cone ! Extension of Type_parameter_pfof for the cone version
     Logical(kind=4) :: do_read_ramses_part_id
     Integer(kind=4) :: shell_first_id
     Integer(kind=4) :: shell_last_id
  End Type Type_parameter_pfof_cone

  Type :: Type_parameter_conecreator !< Type for conecreator input parameters common to particles and cells versions
     Character(len=H5STRLEN) :: input_path
     Character(len=H5STRLEN) :: simulation_name
     Character(len=H5STRLEN) :: cone_input_file
     Character(len=H5STRLEN) :: info_cone_input_file
     Character(len=H5STRLEN) :: info_ramses_input_file
     Integer(kind=4)         :: nfile
     Integer(kind=4)         :: first_file
     Integer(kind=4)         :: last_file
     Real(kind=8)            :: cone_max_radius
     Real(kind=8)            :: cube_size
  End Type Type_parameter_conecreator

  Type, extends(Type_parameter_conecreator) :: Type_parameter_conecreator_part !< Extension of Type_parameter_conecreator for the particles version
     Logical(kind=4) :: do_read_ramses_part_id
     Logical(kind=4) :: do_read_potential
     Logical(kind=4) :: do_read_gravitational_field
  End type Type_parameter_conecreator_part

  Type, extends(Type_parameter_conecreator) :: Type_parameter_conecreator_grav !< Extension of Type_parameter_conecreator for the cells version
     Integer(kind=4) :: nlevel
     Integer(kind=4) :: level_min
  End type Type_parameter_conecreator_grav

  !! Variables read from the metadata common to all hdf5 output files
  Type :: Type_common_metadata       ! Type for metadata common to all hdf5 output files
     Character(len=16) :: created_by
     Integer(kind=4)   :: svn_version
     Character(len=16) :: simulation_code
     Character(len=16) :: particle_type
     Integer(kind=4)   :: constant_mass
     Character(len=16) :: units
     Integer(kind=8)   :: npart_simulation
  End type Type_common_metadata

  !! Variables read from the parameter file
  integer :: nb_proc_fof                                  !number of processor use for the fof algorithm <=> number of file+1 to treat
  integer :: ierr                                         !file reading control : well read (if ==0) or error (if /= 0), ie the version of ramsesis not fill => ramses 3
  character(len=H5STRLEN) :: infile                            !name of the parameters file in input
  character (len=H5STRLEN):: file_fof_root                     !base of the name of the file to treat
  character (len=7)       :: bstring                           !Linking length string (ex:_b01000 for b=0.1 -if nul b=0.2)
  character (len=H5STRLEN):: file_cosmo                        !File for tabulated cosmology
  real :: box_len_Mpc                                          !length of the universe in Mpc
  Type(Type_common_metadata) :: common_metadata
  Type(Type_info_ramses) :: simu_info
  Type(Type_parameter_pfof_snap) :: param_pfof

  ! Variables used in modfunctions.f90
  integer(kind=PRI), dimension(:), allocatable :: haloIdentity      !halo IDs of all halos over a core
  real(8), dimension(:,:), allocatable :: haloCOMX, haloCOMV        !center-of-mass position and velocity for all halos over a core
  real(8), dimension(:),   allocatable :: haloRmax, haloVmax        !maximum extension and particle speed for all halos over a core
  real(8), dimension(:),   allocatable :: haloDispX, haloDispV      !position and velocity dispersion for all halos over a core
  real(8), dimension(:,:), allocatable :: haloCenter_potential      !location of the particle with minimum potential for all halos over a core
  real(8), dimension(:),   allocatable :: haloEpot, haloEkin        !potential and kinetic energy for halos over a core
  real(8), dimension(:,:), allocatable :: haloLambda_prime          !normalized angular momentum for halos over a simulation (first normalisation)
  real(8), dimension(:,:), allocatable :: haloLambda                !normalized angular momentum for halos over a simulation (second normalisation)
  real(8), dimension(:),   allocatable :: haloCos_alpha             !angle between the angular momentum axis and the le vecteur propre du tenseur d'inertie
  real(8), dimension(:,:), allocatable :: haloEin_Vec_Inert_Tens    !eigen vectors of the inertial tensor for all halos over a core
  real(8), dimension(:,:), allocatable :: haloEin_Val_Inert_Tens    !eigen values of the inertial tensor for all halos over a core
  real(8), dimension(:,:), allocatable :: haloCircular_values       !M_circ, Rho_circ, V_circ, R_circ where V_circ is max, for all halos over a core
  real(8), dimension(:,:), allocatable :: haloDelta_values          !M178, R178 for all halos over a core
  real(8), dimension(:),   allocatable :: haloRadial_profile        !central values of radial bins for profiles
  real(8), dimension(:),   allocatable :: haloRadial_profile178     !central values of R/R178 bins for profiles
  real(8), dimension(:,:), allocatable :: haloVelocity_profile      !circular velocity profile for all halos over a core
  real(8), dimension(:,:), allocatable :: haloDensity_profile       !radial density profile for all halos over a core
  real(8), dimension(:,:), allocatable :: haloDensity_profile178    !radial (r/r178) density profile for all halos over a core
  real(8), dimension(:,:), allocatable :: haloDensity_profile_integrated    !radial integrated density profile for all halos over a core
  real(8), dimension(:,:), allocatable :: haloDensity_profile_integrated178 !radial (r/r178) integrated density profile for all halos over a core

  real(8), dimension(:), allocatable :: HaloRvir                    !virial radius for all halos over a core (given in phys)
  real(8), dimension(:), allocatable :: HaloRvircom                 !virial radius for all halos over a core (given in com)
  real(8), dimension(:), allocatable :: HaloVvir                    !virial speed for all halos over a core
  real(8), dimension(:), allocatable :: HaloMvir                    !renormalized mass for all halos over a core
  integer(8), dimension(:), allocatable :: HaloNpart                !number of particles
  real(8), dimension(:), allocatable :: HaloVolume                  !volume of the halo

  real(8), parameter :: msun=1.98892d33, mpc=3.08568d24             !solar mass and megaparsec in cgs
  real(8) :: Rho_crit0=1.879d-29                                    !present-day critical density (h^2 g/cm^3)
  real(8) :: PI=3.14159265359879                                    !PI
  real(8) :: G=6.6731E-08                                           !G in cgs
  real(8) :: Rho_m                                                  !background matter density at redshift z
  real(8) :: mass_part                                              !mass of a ColdDarkMatter particle
  real(8) :: Delta_X                                                !length of the littlest grid for non recovery of particle
  real(8) :: random                                                 !random number used for randm_index
  real(8) :: Hz                                                     !Hubble expansion at z
  real(8) :: cos_alpha                                              !angle between the angular momentum axe and the le vecteur propre du tenseur d'inertie

  logical,parameter::do_lightcone=.false.                           !false for regular ramses output
  logical,parameter::do_periodic_boundary=.true.                    !true for regular ramses output
  logical,parameter::do_center_mass=.true.                         !compute halo density profile centered on center-of-mass or minimum-of-potential

  integer(8), dimension(:),allocatable :: idh                       !Parallel FoF ID of each halo
  integer(8)::idh2                                                  !Temporary id
  integer::nb_halos2,mass_halo2                                     !Temporary number and mass
  real::xcenter2,ycenter2,zcenter2                                  !Temporary center
  real(8)::conv,convdens                                                     !Conversion factor

  !!!The following parameters are used for computing the halo density profile
  !!! in shells and integrated (both in Msun/h / (kpc^3)
  !!!Also computed is the following:
  !!! Maximum circular velocity: V_circ_max=sqrt(GM(<R)/R)
  !!! Radius at which V_circ_max is reached: R_circ_max
  !!! Mass within R_circ_max: M_circ_max(<R_circ_max)
  !!! Overdensity (w.r.t background matter density) within R_circ_max   
  integer(kind=4) :: innermost_particles
  real(kind=8)::Rho_critZ !critical density at redshift z
  real(kind=8)::Delta_matter !overdensity w.r.t background matter density at z
  real(kind=8)::Delta_critZ !overdensity w.r.t cricital density at z
  real(kind=8), parameter::min_r=1.d0 !in kpc
!  real(kind=8), parameter::delta_lr=0.01d0 !D(ln R) is the log R bin
  real(kind=8), parameter::delta_lr=0.1d0 !D(ln R) is the log R bin
  real(kind=8), parameter::min_r178=1d-2 !in R/R_178 units
!  real(kind=8), parameter::delta_lr178=0.02d0 !D(ln R/R178) is the log R/R178 bin
  real(kind=8), parameter::delta_lr178=0.2d0 !D(ln R/R178) is the log R/R178 bin
  real(kind=8)::lmin_r, lmin_r178, M178, R178, Diff
  real(kind=8)::V_circ_max,R_circ_max,M_circ_max,Rho_circ_max
!  integer, parameter::NBIN1=800,NBIN2=300
  integer, parameter::NBIN1=100,NBIN2=30


  integer::NN                                                       !Number of elements in tabulated cosmology
  integer, parameter :: NTOT=10000000
  real(8), parameter :: CLIGHT=299792.d0                            !light speed in km/s
  real(8)::ATEMP,EHUBBLETEMPA2,DUMMY,RTEMP
  real(8)::AEXP_TEMP(NTOT),EHUBBLE_TEMP(NTOT),RDIST_TEMP(NTOT)
  real(8)::AEXP(NTOT),EHUBBLE(NTOT),RDIST(NTOT)


#ifdef WITHMPI
  Integer :: procID
  Integer :: procNB
  Integer :: mpierr
#endif


  Character(len=400), dimension(:), allocatable :: filelist    ! list of the filenames that we must read
  Integer :: myfilelistsize
  Integer :: filelistsize            ! number of halo files in the list
  Character(len=400) :: filenameh5 ! name of the parameter file: and example of the 2 different kind of files can be found with the source files

  
  logical:: use_halo_finder_center_as_center_of_mass=.true.

End Module modvariables
