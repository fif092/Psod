!==============================================================================
! Project: pFoF
! File: pfof_cone/src/pfof_cone.f90
! Copyright Fabrice Roy and Vincent Bouillot (2011)
! Fabrice.Roy@obspm.fr
! 
! This software is a computer program whose purpose is to detect dark matter
! haloes in cosmological simulation with a parallel Friend of Friend algorithm.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and, more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!==============================================================================

!> @file
!! This file contains the main program of the light cone implementation of pFoF.

!> pFOF_cone is a distributed implementation of the Friends-of-Friends halo finder algorithm designed for lightcones simulated with RAMSES. <br>
!! This parallel Friend of Friend implementation is based on a sequential
!! implementation written by Edouard Audit (CEA). <br>
!! It has been written by Fabrice Roy -- CNRS / LUTH (Observatoire de Paris) and  Vincent Bouillot. <br>
!! mail: fabrice.roy@obspm.fr <br>
!! mail: vincent.bouillot@obspm.fr <br>

Program pfof_cone

  Use mpi
  Use modmpicommons
  Use modhdf5
  Use modio
  Use modvariables
  Use modmpicom
  Use modfofpara

  Implicit None

  Integer(kind=4) :: mpierr

  ! Initialization of MPI
  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)

  If(procID==0) Then
     Call title()
  End If

  ! Initialization of HDF5
  Call hdf5_init()

  Call readparameters()

  Call h5readmap()

  Call h5readshellinfo()

  Call h5readparticles()

  Call computeminmax()

  Call fofparacone()

  Call deallocall()

  ! Finalize HDF5
  Call hdf5_finalize()

  ! Finalize MPI
  Call Mpi_Finalize(mpierr)

End Program pfof_cone
