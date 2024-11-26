Module modmap

  Use modvariables
  
  Implicit None

  Integer(kind=4) :: fofgrid_imin, fofgrid_imax
  Integer(kind=4) :: fofgrid_jmin, fofgrid_jmax
  Integer(kind=4) :: fofgrid_kmin, fofgrid_kmax
  Integer(kind=4), dimension(:), allocatable :: ict

Contains

  !=======================================================================
  !> 
  Subroutine exploremap()

    Implicit None

    Integer(kind=4) :: factor
    Integer(kind=4) :: if
    Integer(kind=4) :: fmax
    
    fmax = minval(nctab) / 2

    ! fof cube >= shell cube
    ! first diag: fof cube = shell cube
    ! then we apply a multiplication factor () to the edge
    ! fof cube = 8* shell cube (edge*2)
    ! fof cube = 64* shell cube (edge*4)
    ! fof cube = 216* shell cube (edge*6)
    ! fof cube = 512* shell cube (edge*8)
    ! fof cube = 1000* shell cube (edge*10) 
    ! and so on (if necessary) up to min(nctab(i)/2)
    
    factor = 1
    Call map(factor)
    
    Do if = 1, fmax
       factor = 2*if
       Call map(factor)
    End Do

    
  End Subroutine exploremap

  !=======================================================================
  !> 
  Subroutine map(factor)

    Use modio

    Implicit None
    
    Integer(kind=4), intent(in) :: factor
    Integer(kind=4), dimension(3) :: vertex
    Integer(kind=4) :: vertex_imin, vertex_imax
    Integer(kind=4) :: vertex_jmin, vertex_jmax
    Integer(kind=4) :: vertex_kmin, vertex_kmax
    Integer(kind=4) :: imin, imax, jmin, jmax, kmin, kmax
    Integer(kind=4) :: nbc_i, nbc_j, nbc_k, nbc_is2, nbc_js2, nbc_ks2
    Integer(kind=4) :: ncubefof, volumefof
    Integer(kind=4) :: ifof, jfof, kfof, ind
    Integer(kind=8) :: npartfof
    Integer(kind=4) :: np
    Integer(kind=4) :: pid
    Integer(kind=4) :: procNB
    Integer(kind=4), dimension(:,:), allocatable :: ictable
    Integer(kind=4) :: ip
    Integer(kind=4) :: fofgrid_xmin, fofgrid_ymin, fofgrid_zmin
    Integer(kind=4), dimension(3) :: dims
    Integer(kind=4), dimension(:,:,:), allocatable :: coord2pid
    Integer(kind=4), dimension(:,:), allocatable :: pid2coord


    If(fullsky) Then
       vertex(1) = nctab(1) / 2
       vertex(2) = nctab(2) / 2
       vertex(3) = nctab(3) / 2


       ! coord. des cubes au coin de celui du sommet du cone dans la grille de 
       ! cubes shell
       vertex_imin = vertex(1) - factor/2+1
       vertex_imax = vertex(1) + factor/2
       vertex_jmin = vertex(2) - factor/2+1
       vertex_jmax = vertex(2) + factor/2
       vertex_kmin = vertex(3) - factor/2+1
       vertex_kmax = vertex(3) + factor/2
       
       imin = -nctab(1)/2
       imax =  nctab(1)/2
       jmin = -nctab(2)/2
       jmax =  nctab(2)/2
       kmin = -nctab(3)/2
       kmax =  nctab(3)/2
    
       ! nb de cubes fof dans chaque direction:
       nbc_is2 = (imax-factor/2)/factor
       if(mod(imax-factor/2, factor) /= 0) nbc_is2 = nbc_is2 + 1
       nbc_js2 = (jmax-factor/2)/factor
       if(mod(jmax-factor/2, factor) /= 0) nbc_js2 = nbc_js2 + 1
       nbc_ks2 = (kmax-factor/2)/factor
       if(mod(kmax-factor/2, factor) /= 0) nbc_ks2 = nbc_ks2 + 1
       nbc_i = 2*nbc_is2+1
       nbc_j = 2*nbc_js2+1
       nbc_k = 2*nbc_ks2+1
       
       ! coordonnees de la grill fof dans la grill shell
       fofgrid_imin = vertex_imin - nbc_is2*factor
       fofgrid_imax = vertex_imax + nbc_is2*factor
       fofgrid_jmin = vertex_jmin - nbc_js2*factor
       fofgrid_jmax = vertex_jmax + nbc_js2*factor
       fofgrid_kmin = vertex_kmin - nbc_ks2*factor
       fofgrid_kmax = vertex_kmax + nbc_ks2*factor

       fofgrid_xmin = (imin - (nbc_is2*factor-vertex_imin+1) )
       fofgrid_ymin = (jmin - (nbc_js2*factor-vertex_jmin+1) )
       fofgrid_zmin = (kmin - (nbc_ks2*factor-vertex_kmin+1) )

    Else
       vertex(1) = 1
       vertex(2) = nctab(2) / 2
       vertex(3) = nctab(3) / 2
       ! coord. des cubes au coin de celui du sommet du cone dans la grille de 
       ! cubes shell
       vertex_imin = 1
       vertex_imax = factor
       vertex_jmin = vertex(2)-factor/2+1
       vertex_jmax = vertex(2)+factor/2
       vertex_kmin = vertex(3)-factor/2+1
       vertex_kmax = vertex(3)+factor/2
       
       imin =  1
       imax =  nctab(1)
       jmin = -nctab(2)/2
       jmax =  nctab(2)/2
       kmin = -nctab(3)/2
       kmax =  nctab(3)/2
    
       ! nb de cubes fof dans chaque direction:
       nbc_i = imax/factor
       if(mod(imax,factor) /= 0) nbc_i = nbc_i + 1
       nbc_js2 = (jmax-factor/2)/factor
       if(mod(jmax-factor/2, factor) /= 0) nbc_js2 = nbc_js2 + 1
       nbc_ks2 = (kmax-factor/2)/factor
       if(mod(kmax-factor/2, factor) /= 0) nbc_ks2 = nbc_ks2 + 1
       nbc_j = 2*nbc_js2+1
       nbc_k = 2*nbc_ks2+1
       
       
       ! coordonnees de la grill fof dans la grill shell
       fofgrid_imin = 1
       fofgrid_imax = nbc_i * factor
       fofgrid_jmin = vertex_jmin - nbc_js2*factor
       fofgrid_jmax = vertex_jmax + nbc_js2*factor
       fofgrid_kmin = vertex_kmin - nbc_ks2*factor
       fofgrid_kmax = vertex_kmax + nbc_ks2*factor
       
       fofgrid_xmin = 0
       fofgrid_ymin = (jmin - (nbc_js2*factor-vertex_jmin+1) )
       fofgrid_zmin = (kmin - (nbc_ks2*factor-vertex_kmin+1) )
    End If

    ncubefof = nbc_i * nbc_j * nbc_k
    Allocate(npartcubefof(ncubefof))
    npartcubefof = 0

    Allocate(ictable(factor**3,ncubefof))
    Allocate(ict(factor**3))
    
    Do ifof = 1, nbc_i
       Do jfof = 1, nbc_j
          Do kfof = 1, nbc_k
             ind = ifof + (jfof-1)*nbc_i + (kfof-1)*nbc_i*nbc_j
             Call getnpfromshell(factor,ifof,jfof,kfof, np)
             npartcubefof(ind) = np
             ictable(:,ind) = ict
          End Do
       End Do
    End Do

    procNB = nbc_i * nbc_j * nbc_k - count(npartcubefof==0)

    Allocate(neigh(6,procNB))
    Allocate(pictable(factor**3,procNB))
    Allocate(my_npart_tab(procNB))
    Allocate(boundaries(6,procNB))

    Allocate(coord2pid(nbc_i, nbc_j, nbc_k))
    Allocate(pid2coord(3,procNB))

    coord2pid = 0
    pid2coord = 0
    pid = 1
    !! On ecrit les cartes pour pfof_cone
    Do ifof = 1, nbc_i
       Do jfof = 1, nbc_j
          Do kfof = 1, nbc_k
             ind = ifof + (jfof-1)*nbc_i + (kfof-1)*nbc_i*nbc_j
             If(npartcubefof(ind) /= 0) Then
                ! on assigne un process, 
                coord2pid(ifof,jfof,kfof) = pid
                pid2coord(1,pid) = ifof
                pid2coord(2,pid) = jfof
                pid2coord(3,pid) = kfof
                my_npart_tab(pid) = npartcubefof(ind)
                pictable(:,pid) = ictable(:,ind)
                boundaries(1,pid) = (fofgrid_xmin + (ifof-1)*factor)*cubesize
                boundaries(2,pid) = (fofgrid_xmin + ifof*factor)*cubesize
                boundaries(3,pid) = (fofgrid_ymin + (jfof-1)*factor)*cubesize
                boundaries(4,pid) = (fofgrid_ymin + jfof*factor)*cubesize
                boundaries(5,pid) = (fofgrid_zmin + (kfof-1)*factor)*cubesize
                boundaries(6,pid) = (fofgrid_zmin + kfof*factor)*cubesize
                pid = pid+1
             End If
          End Do
       End Do
    End Do

    Do ip = 1, procNB
       ifof = pid2coord(1,ip)
       jfof = pid2coord(2,ip)
       kfof = pid2coord(3,ip)

       If(ifof==1) Then
          neigh(1,ip) = -1
       Else
          ind = (ifof-1) + (jfof-1)*nbc_i + (kfof-1)*nbc_i*nbc_j
          If(npartcubefof(ind) /= 0) Then
             neigh(1,ip) = coord2pid(ifof-1,jfof,kfof) - 1
          Else
             neigh(1,ip) = -1
          End If
       End If
         
       If(ifof==nbc_i) Then
          neigh(2,ip) = -1
       Else
          ind = (ifof+1) + (jfof-1)*nbc_i + (kfof-1)*nbc_i*nbc_j
          If(npartcubefof(ind) /= 0) Then
             neigh(2,ip) = coord2pid(ifof+1,jfof,kfof) - 1
          Else
             neigh(2,ip) = -1
          End If
       End If

       If(jfof==1) Then
          neigh(3,ip) = -1
       Else
          ind = ifof + (jfof-2)*nbc_i + (kfof-1)*nbc_i*nbc_j
          If(npartcubefof(ind) /= 0) Then
             neigh(3,ip) = coord2pid(ifof,jfof-1,kfof) - 1
          Else
             neigh(3,ip) = -1
          End If
       End If

       If(jfof==nbc_j) Then
          neigh(4,ip) = -1
       Else
          ind = ifof + jfof*nbc_i + (kfof-1)*nbc_i*nbc_j
          If(npartcubefof(ind) /= 0) Then
             neigh(4,ip) = coord2pid(ifof,jfof+1,kfof) - 1
          Else
             neigh(4,ip) = -1
          End If
       End If

       If(kfof==1) Then
          neigh(5,ip) = -1
       Else
          ind = ifof + (jfof-1)*nbc_i + (kfof-2)*nbc_i*nbc_j
          If(npartcubefof(ind) /= 0) Then
             neigh(5,ip) = coord2pid(ifof,jfof,kfof-1) - 1
          Else
             neigh(5,ip) = -1
          End If
       End If

       If(kfof==nbc_k) Then
          neigh(6,ip) = -1
       Else
          ind = ifof + (jfof-1)*nbc_i + kfof*nbc_i*nbc_j
          If(npartcubefof(ind) /= 0) Then
             neigh(6,ip) = coord2pid(ifof,jfof,kfof+1) - 1
          Else
             neigh(6,ip) = -1
          End If
       End If


    End Do

    dims(1) = nbc_i
    dims(2) = nbc_j
    dims(3) = nbc_k
    

    ! volume occupe par les cubes fof en unite de volume d'un cube shell
    volumefof = ncubefof * factor**3

    npartfof = sum(npartcubefof)
    If(npartfof /= npart) Print *,'Error: number of particles in fof cubes differ from number of particles read from the shells!'

    Print *,' '
    Print '(A)','Report:'
    Print '(A,I12,A,I12)','Nb of particles in fof cubes/nb of particles in shell files: ', npartfof,' /',npart
    Print '(A,I7)','Nb of fof cubes and processes: ', procNB
    Print '(A,I10)','Max number of particles in a fof cube: ',maxval(npartcubefof)
    Print '(A,I10)','Min number of particles in a fof cube: ',minval(npartcubefof,1,npartcubefof/=0)

    Call h5writemap(factor, procNB, dims)


    Deallocate(neigh)
    Deallocate(ict)
    Deallocate(ictable)
    Deallocate(pictable)
    Deallocate(npartcubefof)
    Deallocate(my_npart_tab)
    Deallocate(boundaries)
    Deallocate(coord2pid)
    Deallocate(pid2coord)
    
  End Subroutine map


  !=======================================================================
  !> 
  Subroutine getnpfromshell(factor,if,jf,kf,np)

    Implicit None

    Integer(kind=4), intent(in) :: factor, if, jf, kf
    Integer(kind=4) :: np

    Integer(kind=4) :: imin, imax, jmin, jmax, kmin, kmax
    Integer(kind=4) :: ic, jc, kc, ind
    
    np = 0

    imin = fofgrid_imin + (if-1)*factor
    imax = imin + factor - 1
    If(imin < 1) imin = 1
    If(imax > nctab(1)) imax = nctab(1)
    jmin = fofgrid_jmin + (jf-1)*factor
    jmax = jmin + factor - 1
    If(jmin < 1) jmin = 1
    If(jmax > nctab(2)) jmax = nctab(2)
    kmin = fofgrid_kmin + (kf-1)*factor
    kmax = kmin + factor - 1
    If(kmin < 1) kmin = 1
    If(kmax > nctab(3)) kmax = nctab(3)

    ict = 0
    ind = 1
    Do ic = imin, imax
       Do jc = jmin, jmax
          Do kc = kmin, kmax
             ict(ind) = ic + (jc-1)*nctab(1) + (kc-1)*nctab(1)*nctab(2)
             np = np + npartcube(ict(ind))
             ind = ind + 1
          End Do
       End Do
    End Do

  End Subroutine getnpfromshell

End Module modmap
