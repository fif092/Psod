!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modparam

  Character(len=200) :: dirname
  Character(len=195) :: filename
  Character(len=400) :: shellname
  Character(len=200), dimension(:), allocatable :: filelist    ! list of the filenames that we must read
  Integer(kind=4) :: filenb
  Integer(kind=4) :: ffile
  Integer(kind=4) :: lfile
  Integer(kind=8) :: npart
  Integer(kind=4) :: nstride
  Logical(kind=4) :: readramsespartid
  Logical(kind=4) :: fullsky
  Real(kind=8) :: cubesize
  Real(kind=8) :: rmax
  Real(kind=8) :: thetay
  Real(kind=8) :: thetaz


  Namelist/infocone/rmax, thetay, thetaz, npart, fullsky
  Namelist/infofile/dirname, filename, filenb, ffile, lfile, nstride, readramsespartid
  Namelist/infooutput/shellname, cubesize

End Module modparam
