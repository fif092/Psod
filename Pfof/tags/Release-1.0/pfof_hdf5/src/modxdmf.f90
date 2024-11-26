!> @file
!! This file contains subroutine used to write xdmf files corresponding to the HDF5 files written by pFOF.
!! The xdmf files can be read by Paraview to visualize the contents of the HDF5 files.

!> This module contains subroutine used to write xdmf files corresponding to the HDF5 files written by pFOF.
!! The xdmf files can be read by Paraview to visualize the contents of the HDF5 files.
Module modxdmf
  
  Implicit None

Contains
  

  Subroutine writehaloxdmf
    
    

  End Subroutine writehaloxdmf

  !=======================================================================
  !> Writes a single xdmf file containing the references to each HDF5 halo file written by pFOF.
  Subroutine writesinglehaloxdmf(procNB,filebase,halonbtab,halomasstab,haloidtab)
    
    Use modconst

    Integer(kind=4), intent(in) :: procNB
    Character(len=76), intent(in) :: filebase
    Integer(kind=4), dimension(procNB), intent(in) :: halonbtab
    Integer(kind=4), dimension(*), intent(in) :: halomasstab
    Integer(kind=PRI), dimension(*), intent(in) :: haloidtab
    Character(len=85) :: filehdf5
    Character(len=81) :: filexdmf
    Character(len=12) :: npartchar
    Character(len=16) :: groupname
    Character(len=5) :: pchar
    Integer :: p, h, ind

    filexdmf = trim(filebase)//'.xdmf'

    groupname = "halo_00000000000"
    
    Open(Unit=30,File=filexdmf)
    Write(30,'(A)') '<?xml version="1.0" ?>'
    Write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    Write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">'
    Write(30,'(A)') '  <Domain>'
    Write(30,'(A)') '    <Grid Name="Halos collection" GridType="Collection" CollectionType="Spatial">'
    ind = 1
    Do p = 1, procNB
       Write(pchar(1:5),'(I5.5)') p-1
       filehdf5 = trim(filebase)//'_'//pchar//'.h5'
       Write(30,'(A)') '    <Grid Name="'//trim(filehdf5)//'" GridType="Collection" CollectionType="Spatial">'
       Do h = 1, halonbtab(p)
          Write(groupname(6:16),'(I11.11)') haloidtab(ind)
          Write(npartchar(1:12),'(I12)') halomasstab(ind)
          ind = ind + 1
          Write(30,'(A)') '    <Grid Name="'//groupname//'" GridType="Uniform">'
          Write(30,'(A)') '      <Topology TopologyType="Polyvertex" NodesPerElement="'//trim(adjustl(npartchar))//'">'
          Write(30,'(A)') '      </Topology>'
          Write(30,'(A)') '      <Geometry GeometryType="XYZ">'
          Write(30,'(A)') '	<DataItem DataType="Float" Dimensions="'//trim(adjustl(npartchar))//' 3" Format="HDF">'
          Write(30,'(A)') '	  '//trim(filehdf5)//':/'//trim(groupname)//'/pos'
          Write(30,'(A)') '	</DataItem>'
          Write(30,'(A)') '      </Geometry>'
          Write(30,'(A)') '      <Attribute AttributeType="Scalar" Center="Node" Name="ID">'
          Write(30,'(A)') '	<DataItem DataType="Int" Precision="8" Dimensions="'//trim(adjustl(npartchar))//' 1" Format="HDF">'
          Write(30,'(A)') '	  '//trim(filehdf5)//':/'//trim(groupname)//'/ID'
          Write(30,'(A)') '	</DataItem>'
          Write(30,'(A)') '      </Attribute>'
          Write(30,'(A)') '      <Attribute AttributeType="Vector" Center="Node" Name="velocity">'
          Write(30,'(A)') '	<DataItem DataType="Float" Dimensions="'//trim(adjustl(npartchar))//' 3" Format="HDF">'
          Write(30,'(A)') '	  '//trim(filehdf5)//':/'//trim(groupname)//'/vel'
          Write(30,'(A)') '	</DataItem>'
          Write(30,'(A)') '      </Attribute>'
          Write(30,'(A)') '    </Grid>'
       End Do
       Write(30,'(A)') '    </Grid>'
    End Do
    Write(30,'(A)') '    </Grid>'
    Write(30,'(A)') '  </Domain>'
    Write(30,'(A)') '</Xdmf>'
    Close(30)



  End Subroutine writesinglehaloxdmf

  
  !=======================================================================
  !> Writes one xmdf file for each HDF5 cube file.
  Subroutine writecubexdmf(filebase,npart)

    Character(len=82), intent(in) :: filebase
    Integer(kind=4), intent(in) :: npart
    Character(len=85) :: filehdf5
    Character(len=87) :: filexdmf
    Character(len=12) :: npartchar

    Write(npartchar(1:12),'(I12)') npart
    filexdmf = trim(filebase)//'.xdmf'
    filehdf5 = trim(filebase)//'.h5'

    Open(Unit=30,File=filexdmf)
    Write(30,'(A)') '<?xml version="1.0" ?>'
    Write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    Write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">'
    Write(30,'(A)') '  <Domain>'
    Write(30,'(A)') '    <Grid Name="Cube" GridType="Uniform">'
    Write(30,'(A)') '      <Topology TopologyType="Polyvertex" NodesPerElement="'//trim(adjustl(npartchar))//'">'
    Write(30,'(A)') '      </Topology>'
    Write(30,'(A)') '      <Geometry GeometryType="XYZ">'
    Write(30,'(A)') '	<DataItem DataType="Float" Dimensions="'//trim(adjustl(npartchar))//' 3" Format="HDF">'
    Write(30,'(A)') '	  '//trim(filehdf5)//':/pos'
    Write(30,'(A)') '	</DataItem>'
    Write(30,'(A)') '      </Geometry>'
    Write(30,'(A)') '      <Attribute AttributeType="Scalar" Center="Node" Name="ID">'
    Write(30,'(A)') '	<DataItem DataType="Int" Precision="8" Dimensions="'//trim(adjustl(npartchar))//' 1" Format="HDF">'
    Write(30,'(A)') '	  '//trim(filehdf5)//':/ID'
    Write(30,'(A)') '	</DataItem>'
    Write(30,'(A)') '      </Attribute>'
    Write(30,'(A)') '      <Attribute AttributeType="Vector" Center="Node" Name="velocity">'
    Write(30,'(A)') '	<DataItem DataType="Float" Dimensions="'//trim(adjustl(npartchar))//' 3" Format="HDF">'
    Write(30,'(A)') '	  '//trim(filehdf5)//':/vel'
    Write(30,'(A)') '	</DataItem>'
    Write(30,'(A)') '      </Attribute>'
    Write(30,'(A)') '    </Grid>'
    Write(30,'(A)') '  </Domain>'
    Write(30,'(A)') '</Xdmf>'
    Close(30)


  End Subroutine writecubexdmf

  !=======================================================================  
  !> Writes a single xdmf file containing the references to each HDF5 cube file written by pFOF.
  Subroutine writesinglecubexdmf(procNB,filebase,nparttab)
    
    Integer(kind=4), intent(in) :: procNB
    Character(len=76), intent(in) :: filebase
    Integer(kind=4), dimension(procNB), intent(in) :: nparttab
    Character(len=85) :: filehdf5
    Character(len=81) :: filexdmf
    Character(len=12) :: npartchar
    Character(len=5) :: pchar
    Integer :: p

    filexdmf = trim(filebase)//'.xdmf'

    Open(Unit=30,File=filexdmf)
    Write(30,'(A)') '<?xml version="1.0" ?>'
    Write(30,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    Write(30,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">'
    Write(30,'(A)') '  <Domain>'
    Write(30,'(A)') '    <Grid Name="Cube collection" GridType="Collection" CollectionType="Spatial">'
    Do p = 1, procNB
       Write(pchar(1:5),'(I5.5)') p-1
       Write(npartchar(1:12),'(I12)') nparttab(p)
       filehdf5 = trim(filebase)//'_'//pchar//'.h5'
       Write(30,'(A)') '    <Grid Name="Cube'//pchar//'" GridType="Uniform">'
       Write(30,'(A)') '      <Topology TopologyType="Polyvertex" NodesPerElement="'//trim(adjustl(npartchar))//'">'
       Write(30,'(A)') '      </Topology>'
       Write(30,'(A)') '      <Geometry GeometryType="XYZ">'
       Write(30,'(A)') '	<DataItem DataType="Float" Dimensions="'//trim(adjustl(npartchar))//' 3" Format="HDF">'
       Write(30,'(A)') '	  '//trim(filehdf5)//':/pos'
       Write(30,'(A)') '	</DataItem>'
       Write(30,'(A)') '      </Geometry>'
       Write(30,'(A)') '      <Attribute AttributeType="Scalar" Center="Node" Name="ID">'
       Write(30,'(A)') '	<DataItem DataType="Int" Precision="8" Dimensions="'//trim(adjustl(npartchar))//' 1" Format="HDF">'
       Write(30,'(A)') '	  '//trim(filehdf5)//':/ID'
       Write(30,'(A)') '	</DataItem>'
       Write(30,'(A)') '      </Attribute>'
       Write(30,'(A)') '      <Attribute AttributeType="Vector" Center="Node" Name="velocity">'
       Write(30,'(A)') '	<DataItem DataType="Float" Dimensions="'//trim(adjustl(npartchar))//' 3" Format="HDF">'
       Write(30,'(A)') '	  '//trim(filehdf5)//':/vel'
       Write(30,'(A)') '	</DataItem>'
       Write(30,'(A)') '      </Attribute>'
       Write(30,'(A)') '    </Grid>'
    End Do
    Write(30,'(A)') '    </Grid>'
    Write(30,'(A)') '  </Domain>'
    Write(30,'(A)') '</Xdmf>'
    Close(30)


  End Subroutine writesinglecubexdmf


End Module modxdmf
