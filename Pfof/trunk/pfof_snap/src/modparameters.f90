! Copyright (c) 2007-2016 CNRS, Fabrice Roy
! Author: Fabrice Roy (LUTH/CNRS/PSL), fabrice.roy@obspm.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

!> @file
!! Input parameters for pfof_snap.
!! @brief
!! The parameters are stored in a user-defined data type Type_parameter_pfof_snap.<br>
!! This datatype is defined in common/src/modconstant.f90
!> @author Fabrice Roy 

!> Input parameters for pfof_snap. <br>
Module modparameters

  Use modconstant
  Implicit None

  Type(Type_parameter_pfof_snap) :: param !< Input parameters for pfof_snap.

End Module modparameters
