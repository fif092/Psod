!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

!=======================================================================
! Generic program reading a set of hdf5 halo files produced by pFoF
! and analyzing them
! You can analyze either a subset of halos by providing thier ID (analyzelist)
! or each halo in a subset of files by providing the name of the files (analyzefile)
!=======================================================================
Program haloanalyzer

  Use modanalyzeengines
  Use modfunctions
  Use modiocommons
  Use modvariables
  use modreadhalo

  Implicit none

  Type(ListOfVar)    :: var      ! list of the datasets needed for your analyzis
  Character(len=400) :: filename ! name of the parameter file: and example of the 2 different kind of files can be found with the source files

  Integer :: narg, iargc
  Integer :: haloNB_glob
  Integer(kind=8) :: npart8
  Real(8) :: mass_CDM_particle    !mass of a ColdDarkMatter particle
  Character(len=H5STRLEN) :: adata
  Character(len=50) :: origin
  Character(len=15) :: codeversion
  
  Character(len=400) :: hf_filenamenml,hf_filename !< Name of the file
  Namelist / File / hf_filename


#ifdef WITHMPI
  Call Mpi_Init(mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,procNB,mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,procID,mpierr)
#endif

  Call hdf5_init()
  
  narg = iargc()
  if(narg .lt. 1)then
    if(procID==0) then
       print*,'You must use this program with a parameter file in argument'
       print*,'If you need help, type "help" in argument'
    endif
    call MPI_FINALIZE(mpierr)
    stop
  end if
  call getarg(1,infile)

  if (trim(infile) == 'help')then
     if(procID==0) then
        !call help
     endif
     call MPI_FINALIZE(mpierr)
     stop
  endif
        
  !You can define the datasets you need to read to perform your analyzis
  var%pos = .true.    ! default is .true.
  var%vel = .true.    ! default is .true.
  var%pot = .false.    ! default is .true.
  var%id  = .true.    ! default is .true.

  if(do_center_mass)then
     if(procID==0) then
        print*,'remark: you are using do_center_mass=.true. for profile this is not the default version'
     endif
  endif

  !If use_halo_finder_center_as_center_of_mass=.true., force to use halo_finder_center. Filename hardcoded, slow, memory consuming, and inefficient I/O version -> to be improved. 
  !Just used for testing at the moment. Default should be .false.
  if(use_halo_finder_center_as_center_of_mass)then
     if(procID==0) then
        print*,'WARNING  USE_HALO_FINDER_CENTER_AS_CENTER_OF_MASS  IS TRUE: BETA VERSION'
        print*,'Will replace center of mass by halo finder center read in hfprop file'
     endif
!     hf_filename="/data/yrasera/prepa4096/boxlen82.03125_n128_lcdmw7/post/sod/test_00011/psod_halo_snap_part_hfprop_boxlen82.03125_n128_lcdmw7.h5" !to improve
     hf_filenamenml="fileprop.nml"

     Open(Unit=11,file=trim(hf_filenamenml),status='old')
     Read(11,nml=File)
     Close(11)
     if(procID==0) then
        print*,'filenameprop',hf_filename
     endif
     call h5readhalomass(hf_filename)
     !print*,'haloID',minval(hf_haloID),maxval(hf_haloID)
     !print*,'mass',minval(hf_halomass),maxval(hf_halomass)
     !print*,'pos',minval(hf_halocompos),maxval(hf_halocompos)
     !print*,'vel',minval(hf_halocomvel),maxval(hf_halocomvel)
  endif




  !The first argument is the name of the file containing the list of the ID of the halos that you want to analyze
  !The second argument is the name of your analyzis subroutine (defined in modfunctions.f90) 
  !compos is an example that you can find in modfunctions.f90
  !The third argument (optional, used in the following example) is the list of variables (datasets) that you need to read from the HDF5 files
  !filename = 'halolist.nml'
  !Call analyzelist(filename,compos,var)
  
  !The first argument is the name of the file containing the list of the files that you want to analyze
  !The second argument is the name of your analyzis subroutine (defined in modfunctions.f90) 
  !comvel is an example that you can find in modfunctions.f90
  !The third argument (optional, not used in the following example) is the list of variables (datasets) that you need to read from the HDF5 files
  filename = 'filelist.nml'
  Call analyzefile(filename, halo_properties, var)

  Call Mpi_AllReduce(nbhaloanalyzed,haloNB_glob,1,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

  if(use_halo_finder_center_as_center_of_mass)then
     deAllocate(hf_halocomPos, hf_halocomVel, hf_haloID, hf_haloMass)
  endif

  !Example of parallel hdf5 output
  !Use parallel HDF5 only if you compile with -DWITHHDF5
  Open(Unit=10, File='haloanalyzer.version', status='old', Iostat=ierr)
  If(ierr/=0) Then
     print*, 'version of haloanalyzer not defined: file not found'
     codeversion = 'undefined'
  Else
     Read(10,*,Iostat=ierr) codeversion
     If(ierr/=0) Then
        print*, 'version of haloanalyzer not defined: file empty'
        codeversion='undefined'
     End If
     Close(10)
  End If
  origin = 'Created with haloanalyzer version '//codeversion

#ifdef SOD
  filenameh5 = 'psod_halo_snap_part_prop_'//trim(param_pfof%simulation_name)//'.h5'
#else
  filenameh5 = 'pfof_halo_snap_part_prop_'//trim(param_pfof%simulation_name)//'.h5'
#endif

 Call hdf5_create_mpi_file(filenameh5, Mpi_Comm_World, file_id)
 name = 'haloNB'
 Call hdf5_write_attr(file_id, name, haloNB_glob)

 groupname = 'data'
 Call hdf5_create_group(file_id, groupname, gr_id)
 
groupname = 'properties'
 Call hdf5_create_group(gr_id, groupname, prop_id)

 name = 'identity_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloIdentity), haloIdentity, Mpi_Comm_World)
 name = 'number_particles_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloNpart), haloNpart, Mpi_Comm_World)
 name = 'mass_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloMvir),   haloMvir, Mpi_Comm_World)
 name = 'normalisation_mass_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloMvir),   haloMvir,  Mpi_Comm_World)
 name = 'normalisation_radius_physical_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloRvir),   haloRvir,  Mpi_Comm_World)
 name = 'normalisation_radius_comoving_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloRvircom),   haloRvircom,  Mpi_Comm_World)
 name = 'normalisation_velocity_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloVvir),   haloVvir,  Mpi_Comm_World)
 name = 'position_center_of_mass_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloCOMX,2), haloCOMX,  Mpi_Comm_World)
 name = 'velocity_center_of_mass_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloCOMV,2), haloCOMV,  Mpi_Comm_World)
 name = 'radius_maximum_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloRmax),   haloRmax,  Mpi_Comm_World)
 name = 'velocity_maximum_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloVmax),   haloVmax,  Mpi_Comm_World)
 name = 'dispersion_position_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloDispX),  haloDispX, Mpi_Comm_World)
 name = 'dispersion_velocity_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloDispV),  haloDispV, Mpi_Comm_World)
 name = 'energy_self_binding_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloEpot),   haloEpot,  Mpi_Comm_World)
 name = 'energy_kinetic_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloEkin),   haloEkin,  Mpi_Comm_World)
 name = 'angular_momentum_lambda_prime_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloLambda_prime,2), haloLambda_prime, Mpi_Comm_World)
 name = 'angular_momentum_lambda_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloLambda,2), haloLambda, Mpi_Comm_World)
 name = 'inertia_minor_axis_vector_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloEin_Vec_Inert_Tens,2), haloEin_Vec_Inert_Tens, Mpi_Comm_World)
 name = 'inertia_eigen_values_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloEin_Val_Inert_Tens,2), haloEin_Val_Inert_Tens, Mpi_Comm_World)
 name = 'cosine_angle_angular_momentum_minor_axis_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloCos_alpha), haloCos_alpha, Mpi_Comm_World)
 name = 'position_most_bounded_halo'
 Call hdf5_write_mpi_data(prop_id, name, 3, size(haloCenter_potential,2), haloCenter_potential, Mpi_Comm_World)
 name = 'circular_values_halo'
 Call hdf5_write_mpi_data(prop_id, name, 4, size(haloCircular_values,2), haloCircular_values, Mpi_Comm_World)
 name = 'delta_values_halo'
 Call hdf5_write_mpi_data(prop_id, name, 2, size(haloDelta_values,2), haloDelta_values, Mpi_Comm_World)
 name = 'volume_halo'
 Call hdf5_write_mpi_data(prop_id, name,    size(haloVolume),   haloVolume,  Mpi_Comm_World)

 Call hdf5_close_group(prop_id)

 groupname = 'profiles'
 Call hdf5_create_group(gr_id, groupname, prof_id)

 name = 'profile_radial_bins_halo'
 Call hdf5_write_data(prof_id, name, NBIN1, haloRadial_profile)
 name = 'profile_radial_bins_rdelta_halo'
 Call hdf5_write_data(prof_id, name, NBIN2, haloRadial_profile178)
 name = 'profile_circular_velocity_halo'
 Call hdf5_write_mpi_data(prof_id, name, NBIN1, size(haloVelocity_profile,2), haloVelocity_profile, Mpi_Comm_World)
 name = 'profile_density_halo'
 Call hdf5_write_mpi_data(prof_id, name, NBIN1, size(haloDensity_profile,2), haloDensity_profile, Mpi_Comm_World)
 name = 'profile_density_integrated_halo'
 Call hdf5_write_mpi_data(prof_id, name, NBIN1, size(haloDensity_profile_integrated,2), haloDensity_profile_integrated,Mpi_Comm_World)
 name = 'profile_density_bins_rdelta_halo'
 Call hdf5_write_mpi_data(prof_id, name, NBIN2, size(haloDensity_profile178,2), haloDensity_profile178, Mpi_Comm_World)
 name = 'profile_density_integrated_bins_rdelta_halo'
 Call hdf5_write_mpi_data(prof_id,name, NBIN2,size(haloDensity_profile_integrated178,2),&
                        & haloDensity_profile_integrated178,Mpi_Comm_World)

 Call hdf5_close_group(prof_id)
 Call hdf5_close_group(gr_id)

    npart8 = 2**(simu_info%levelmin)
    npart8 = npart8**3

    Call h5write_meta_common(file_id, npart8)
    Call h5write_meta_pfof_parameter(file_id, param_pfof)
    Call h5write_meta_info_ramses(file_id, simu_info)
    Call h5write_info_haloanalyzer(file_id)

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)
    name = 'halo_finder'
#ifdef SOD
    adata = 'psod'
#else
    adata = 'pfof'
#endif
    Call hdf5_write_attr(gr_id, name, adata)
    ! Write type as attribute    
    name = 'file_type'
    adata = 'halo_properties'
    Call hdf5_write_attr(gr_id, name, adata)
    ! write the number of halos as an attribute
    name = 'nhalo_file'
    Call hdf5_write_attr(gr_id, name, haloNB_glob)
    ! Write do_center_mass as an integer attribute (1=true, 0=false)
    name = 'center_type'
    if(do_center_mass)then
       adata = 'center_of_mass'
    else
       adata = 'minimum_of_potential'
    endif
    Call hdf5_write_attr(gr_id, name, adata)
    ! Write mass of a ColdDarkMatter particle as an attribute
    mass_CDM_particle = mass_part*simu_info%unit_d*simu_info%unit_l**3          !in grams
    mass_CDM_particle = real(mass_CDM_particle/msun * simu_info%h0/100.,kind=8) !in Msun/h
    name = 'mass_CDM_particle'
    Call hdf5_write_attr(gr_id, name, mass_CDM_particle)
    Call hdf5_close_group(gr_id)

 Call hdf5_close_mpi_file(file_id)

  If(Allocated(haloIdentity))           Deallocate(haloIdentity)
  If(Allocated(haloCOMX))               Deallocate(haloCOMX)
  If(Allocated(haloCOMV))               Deallocate(haloCOMV)
  If(Allocated(haloRmax))               Deallocate(haloRmax)
  If(Allocated(haloVmax))               Deallocate(haloVmax)
  If(Allocated(haloDispX))              Deallocate(haloDispX)
  If(Allocated(haloDispV))              Deallocate(haloDispV)
  If(Allocated(haloMvir))               Deallocate(haloMvir)
  If(Allocated(haloNpart))              Deallocate(haloNpart)
  If(Allocated(haloRvir))               Deallocate(haloRvir)
  If(Allocated(haloRvircom))            Deallocate(haloRvircom)
  If(Allocated(haloVvir))               Deallocate(haloVvir)
  If(Allocated(haloEpot))               Deallocate(haloEpot)
  If(Allocated(haloEkin))               Deallocate(haloEkin)
  If(Allocated(haloLambda_prime))       Deallocate(haloLambda_prime)
  If(Allocated(haloLambda))             Deallocate(haloLambda)
  If(Allocated(haloEin_Vec_Inert_Tens)) Deallocate(haloEin_Vec_Inert_Tens)
  If(Allocated(haloEin_Val_Inert_Tens)) Deallocate(haloEin_Val_Inert_Tens)
  If(Allocated(haloCos_alpha))          Deallocate(haloCos_alpha)
  If(Allocated(haloCenter_potential))   Deallocate(haloCenter_potential)
  If(Allocated(haloCircular_values))    Deallocate(haloCircular_values)
  If(Allocated(haloDelta_values))       Deallocate(haloDelta_values)
  If(Allocated(haloRadial_profile))     Deallocate(haloRadial_profile)
  If(Allocated(haloRadial_profile178))  Deallocate(haloRadial_profile178)
  If(Allocated(haloVelocity_profile))   Deallocate(haloVelocity_profile)
  If(Allocated(haloDensity_profile))    Deallocate(haloDensity_profile)
  If(Allocated(haloDensity_profile178)) Deallocate(haloDensity_profile178)
  If(Allocated(haloDensity_profile_integrated))    Deallocate(haloDensity_profile_integrated)
  If(Allocated(haloDensity_profile_integrated178)) Deallocate(haloDensity_profile_integrated178)
  If(Allocated(haloVolume))               Deallocate(haloVolume)
  Call hdf5_finalize()

#ifdef WITHMPI
  Call Mpi_Barrier(MPI_COMM_WORLD,mpierr)
  if(procID==0) print*,'Run completed'
  Call Mpi_Finalize(mpierr)
#endif

End Program haloanalyzer
