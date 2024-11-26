!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================
Module modfunctions
  
  Use modvariable

  !! WARNING: you can declare here the output variables you may want to deal with in the 
  !! main program, for instance if you want to write them in an output file.
  !! If the variable is an allocatable array, you have to allocate it the first
  !! time you use your analyze function (see compos for an example).
  !! But you have to deallocate these arrays before a 2nd call to analyze_list or analyze_all 
  !! if you want to call 2 different analyzes in the same main program.

  Real(kind=8), dimension(:,:), allocatable :: cmpos
  Real(kind=8), dimension(:,:), allocatable :: cmvel
  

Contains

  Subroutine compos() 
    
    Integer(kind=4) :: ip
    Integer(kind=4) :: np
    Real(kind=8), dimension(3) :: x
    
    If(.not.Allocated(cmpos)) Then 
       Allocate(cmpos(3,nbhaloanalyzed))
    End If

    np = size(pos,2)
    
    x = 0.d0
    Do ip = 1, np
       x = x + pos(:,ip)
    End Do
    
    cmpos(:,currenthalo) = x / real(np)
    
  End Subroutine  compos



  !=======================================================================
  
  Subroutine comvel() 
    
    Integer(kind=4) :: ip
    Integer(kind=4) :: np
    Real(kind=8), dimension(3) :: v
    
    If(.not.Allocated(cmvel)) Allocate(cmvel(3,nbhaloanalyzed))

    np = size(vel,2)
    
    v = 0.d0
    Do ip = 1, np
       v = v + vel(:,ip)
    End Do
    cmvel(:,currenthalo) = v / real(np)
    
  End Subroutine  comvel
  
  
End Module modfunctions
