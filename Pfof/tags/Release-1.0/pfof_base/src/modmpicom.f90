Module modmpicom

  Use modconst

  Character, dimension(:),allocatable :: header
  Integer(kind=4) :: mpierr
  Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqr1,mpireqr2,mpireqr3,mpireqr4
  Integer(kind=4) :: procID, procNB
  Integer(kind=4) :: h_pos, h_length
  Integer(kind=4) :: dims(3)
  Integer(kind=4) :: MPICube
  Integer(kind=4) :: CubeCoord(3)
  Integer(kind=4) :: voisin(6)
  Logical :: periods(3)
  
Contains

  Subroutine EmergencyStop(message,errcode)

    Character(len=*), Intent(in) :: message
    Integer, Intent(in) :: errcode

    Write(*,1000) message,procID

1000 Format('*** ',A,' on process ',I5.5,'. ***')

    If(procID==0) Close(50)
    Close(54)

    Call Mpi_Abort(MPICube,errcode,mpierr)
    Stop

  End Subroutine EmergencyStop

End Module modmpicom
