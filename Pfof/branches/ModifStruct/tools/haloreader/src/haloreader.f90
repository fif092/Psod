Program haloreader

  Use modconstant
  Use modreadhalo
  Use modhdf5

  Implicit none

  Real(kind=8), dimension(:,:), allocatable :: position, velocity
  Real(kind=8), dimension(:), allocatable :: rmax
  Integer(kind=PRI), dimension(:), allocatable :: id
  Integer(kind=4), dimension(:), allocatable :: np
  Type(Type_parameter_pfof_cone) :: pfof
  Type(Type_info_ramses) :: ramses
  Type(Type_info_cone_part) :: cone
  Type(Type_common_metadata) :: common_metadata
  Character(len=400) :: filename 
  
  Call hdf5_init()

  filename = '/home/roy/Travail/Devel/Cosmologie/test/trunk/pfof_halo_cone_part_hfprop_testtrunk.h5'

  Call h5read_halo_hfprop(filename, common_metadata, position_halo=position, velocity_halo=velocity,&
       parameter_pfof=pfof,info_ramses=ramses, info_cone=cone, rmax_halo=rmax, &
       identity_halo=id, npart_halo=np) !  

#ifdef DEBUG
  Print *,'Controle:'
  Print *,'common meta: ',common_metadata
  Print *,'1st position: ',position(:,1)
  Print *,'last position: ',position(:,size(position,2))
  Print *,'1st velocity: ',velocity(:,1)
  Print *,'last velocity:', velocity(:,size(velocity,2))
  Print *,'pfof parameter: ', pfof
#endif

  Call hdf5_finalize()

End Program haloreader
