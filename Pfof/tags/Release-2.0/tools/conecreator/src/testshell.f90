Program testshell

  Use modparam
  Use modio
  Use modvariable
  Use mpi
  Use modsortpart
  Use modhdf5

  Implicit None

  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)

  Call hdf5_init()

  Call readparameters

  Call readramsesfiles

  Call writeramsestext

  Call readwritehdf5

!  Call writehdf5text

  Deallocate(pos, vel, id)
  Deallocate(nparttab)

  Call hdf5_finalize()
  Call Mpi_Finalize(mpierr)

End Program testshell

Subroutine writeramsestext
  
  Use modio
  Use modvariable
  Implicit None

  Character(len=2) :: charpid
  Character(len=12) :: filetext
  Integer :: ip

  Write(charpid(1:2),'(I2.2)') procID

  filetext='ramses'//charpid//'.txt'
  Open(Unit=10, File=filetext)
  Do ip = 1, npartloc
     Write(10,'(6E16.8)')  pos(1,ip), pos(2,ip), pos(3,ip), vel(1,ip), vel(2,ip), vel(3,ip)
  End Do
  Close(10)

End Subroutine writeramsestext


Subroutine readwritehdf5
  
  Use modhdf5
  Use modvariable
  Use modparam

  Implicit None

  Character(len=16) :: dataname, groupname
  Character(len=8) :: charic, filetext
  Integer :: n, ic, ip, deb, fin, npartloc
  Integer, dimension(3) :: nctab

  Integer(kind=hid_t) :: file_id, gr_id

  Call hdf5_open_file(shellname, file_id)
  
  dataname = 'nctab'
  n = 3
  Call hdf5_read_attr(file_id, dataname, n, nctab)

  dataname = 'npartcube'
  n = nctab(1)*nctab(2)*nctab(3)
  If(Allocated(npartcube) ) Then
     Deallocate(npartcube)
  End If
  Allocate(npartcube(n))
  Call hdf5_read_data(file_id, dataname, n, npartcube)

  npartloc = sum(npartcube)
  If(Allocated(pos)) Then
     Deallocate(pos)
  End If
  If(Allocated(vel)) Then
     Deallocate(vel)
  End If
  Allocate(pos(3,npartloc), vel(3, npartloc))

  deb = 1
  Do ic = 1, n
     If(npartcube(ic) /= 0) Then
        Write(charic(1:8),'(I8.8)') ic
        groupname = 'cube'//charic
        Call hdf5_open_group(file_id, groupname, gr_id)

        dataname = 'pos'
        n = npartcube(ic)
        fin = deb + n - 1

        Call hdf5_read_data(gr_id, dataname, 3, n, pos(:,deb:fin))
        dataname = 'vel'
        Call hdf5_read_data(gr_id, dataname, 3, n, vel(:,deb:fin))
        deb = fin + 1

        Call hdf5_close_group(gr_id)
     End If
  End Do

  filetext='hdf5.txt'
  Open(Unit=10, File=filetext)
  Do ip = 1, npartloc
     Write(10,'(6E16.8)')  pos(1,ip), pos(2,ip), pos(3,ip), vel(1,ip), vel(2,ip), vel(3,ip)
  End Do
  Close(10)
  
End Subroutine readwritehdf5
