Module modfunctions
  
  Use modvariables

  !! WARNING: you can declare here the output variables you may want to deal with in the 
  !! main program, for instance if you want to write them in an output file.
  !! If the variable is an allocatable array, you have to allocate it the first
  !! time you use your analyze function (see compos for an example).
  !! But you have to deallocate these arrays before a 2nd call to analyze_list or analyze_all 
  !! if you want to call 2 different analyzes in the same main program.

Contains

  FUNCTION ZDIST(DIST) result(ZD)
  IMPLICIT NONE
  INTEGER i
  REAL*8 DIST,ADIST,ZD

  i=3 !!!!!!WARNING: ramses_input.dat is not accurate enough very close to a=1 so I start at i=3 instead of 2!!!
  do while(RDIST(i)<DIST .and. i<NN)
     i=i+1
  end do
  !Interpolate
  ADIST = AEXP(i)*(DIST-RDIST(i-1))/(RDIST(i)-RDIST(i-1))+ &
        & AEXP(i-1)*(DIST-RDIST(i))/(RDIST(i-1)-RDIST(i))

  !ADIST=AEXP(J-1)+(AEXP(J)-AEXP(J-1))*(DIST-RDIST(J-1))/ (RDIST(J)-RDIST(J-1))
  ZD = 1./ADIST -1
END FUNCTION ZDIST
  !=======================================================================

  FUNCTION E_HUBBLE(ATEMP) result(EH)
  IMPLICIT NONE
  INTEGER i
  REAL*8 ATEMP,EH

  i=3 !!!!!!WARNING: ramses_input.dat is not accurate enough very close to a=1 so I start at i=3 instead of 2!!!
  do while(AEXP(i)>ATEMP .and. i<NN)
     i=i+1
  end do
  !Interpolate
  EH = EHUBBLE(i)*(simu_info%aexp-AEXP(i-1))/(AEXP(i)-AEXP(i-1))+ &
     & EHUBBLE(i-1)*(simu_info%aexp-AEXP(i))/(AEXP(i-1)-AEXP(i))
END FUNCTION E_HUBBLE
  !=======================================================================

  Subroutine halo_properties()
    use modreadhalo !only for use_halo_finder_center_as_center_of_mass mode

    real(kind=8) :: dtemp, zz
    integer(kind=4) :: i,j,np,count_cpt                               !used
    real(kind=8) :: M_vir, R_vir, V_vir                               !virial mass, radius and velocity of the current halo (in cgs)
    real(kind=8), dimension(3) :: comx, comv,hfcomx                          !position and velocity of center-of-mass of halo
    real(kind=8) :: R_max, R_tmp                                      !maximum extension of halo
    real(kind=8) :: V_max, V_tmp                                      !maximum speed of particles within halo
    real(kind=8) :: DispX, DispV                                      !position and velocity dispersion of particles in the current halo
    real(kind=8) :: Epot, Ekin, Etot                                  !potential, kinetic and total energy of the current halo
    real(kind=8) :: Epot_i                                            !potential energy of i-th particle of the current halo
    real(kind=8) :: Epot_max                                          !potential energy of the most tightly bound particle of the current halo
    real(kind=8) :: jtot                                              !total angular momentum of the current halo
    real(kind=8), dimension(3) :: center_pot                          !location of the particle of the curent halo with minimum potential
    real(kind=8), dimension(3) :: ang_momentum                        !angular momentum of the current halo
    real(kind=8), dimension(3) :: Lambda_prime, Lambda                !normalized angular momentum (1st and 2nd normalization)
    real(kind=8), dimension(3,3) :: Inert_Tens, inert_tmp             !inertial tensor of the current halo
    real(kind=8), dimension(3,3) :: Ein_Vec_Inert_Tens                !eigen vectors of the inertial tensor of the current halo
    real(kind=8), dimension(3)   :: Ein_Val_Inert_Tens                !eigen values of the inertial tensor of the current halo
    real(kind=8), dimension(:,:),allocatable :: rand_positions        !array of randomly chosen positions inside a halo
    integer,      dimension(:),  allocatable :: random_index_arr      !array containing all random index for the fast potential energy computation

    integer :: ok                                                     !stat number
    integer :: threshold                                              !limit number of particles for potential energy computation
    integer :: threshold_calc                                         !number of computed particle for potential energy computation
    integer :: random_index                                           !random index for the large halos energy computation
    integer :: nrot                                                   !number of rotation for the search of eigen vector and values
 
    logical :: test

    real(kind=8), dimension(:),allocatable :: density_profile, density_profile178
    real(kind=8), dimension(:),allocatable :: density_profile_integrated, density_profile_integrated178
    real(kind=8), dimension(:),allocatable :: V_circ
    integer, dimension(:),allocatable :: dndlnR, dndlnR178 !D(n)/D(ln R)
    
    logical::ascii=.false.,compute_volume
    integer::jstart
    integer::ihf

    real(kind=8)::bfof

    lmin_r=log(min_r)
    !if(procID==0) print*, lmin_r, NBIN1

    lmin_r178=log(min_r178)
    !if(procID==0) print*, lmin_r178, NBIN2

    !Allocate
    allocate(density_profile(NBIN1),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the density_profile array'
      call MPI_FINALIZE(mpierr)
      stop
    endif
    allocate(density_profile178(NBIN2),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the density_profile array'
      call MPI_FINALIZE(mpierr)
      stop
    endif
    allocate(density_profile_integrated(NBIN1),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the density_profile_integrated array'
      call MPI_FINALIZE(mpierr)
      stop
    endif
    allocate(density_profile_integrated178(NBIN2),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the density_profile_integrated178 array'
      call MPI_FINALIZE(mpierr)
      stop
    endif
    allocate(V_circ(NBIN1),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the V_circ array'
      call MPI_FINALIZE(mpierr)
      stop
    endif
    allocate(dndlnR(NBIN1),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the dndlnR array'
      call MPI_FINALIZE(mpierr)
      stop
    endif
    allocate(dndlnR178(NBIN2),stat=ok)
    if(ok>0)then
      print*,procID, 'Error in allocation of the dndlnR178 array'
      call MPI_FINALIZE(mpierr)
      stop
    endif

    If(.not.Allocated(haloIdentity))           Allocate(haloIdentity(nbhaloanalyzed))
    If(.not.Allocated(haloCOMX))               Allocate(haloCOMX(3,nbhaloanalyzed))
    If(.not.Allocated(haloCOMV))               Allocate(haloCOMV(3,nbhaloanalyzed))
    If(.not.Allocated(haloRmax))               Allocate(haloRmax(nbhaloanalyzed))
    If(.not.Allocated(haloVmax))               Allocate(haloVmax(nbhaloanalyzed))
    If(.not.Allocated(haloDispX))              Allocate(haloDispX(nbhaloanalyzed))
    If(.not.Allocated(haloDispV))              Allocate(haloDispV(nbhaloanalyzed))
    If(.not.Allocated(haloMvir))               Allocate(haloMvir(nbhaloanalyzed))
    If(.not.Allocated(haloRvir))               Allocate(haloRvir(nbhaloanalyzed))
    If(.not.Allocated(haloRvircom))            Allocate(haloRvircom(nbhaloanalyzed))
    If(.not.Allocated(haloVvir))               Allocate(haloVvir(nbhaloanalyzed))
    If(.not.Allocated(haloEpot))               Allocate(haloEpot(nbhaloanalyzed))
    If(.not.Allocated(haloEkin))               Allocate(haloEkin(nbhaloanalyzed))
    If(.not.Allocated(haloLambda_prime))       Allocate(haloLambda_prime(3,nbhaloanalyzed))
    If(.not.Allocated(haloLambda))             Allocate(haloLambda(3,nbhaloanalyzed))
    If(.not.Allocated(haloEin_Vec_Inert_Tens)) Allocate(haloEin_Vec_Inert_Tens(3,nbhaloanalyzed))
    If(.not.Allocated(haloEin_Val_Inert_Tens)) Allocate(haloEin_Val_Inert_Tens(3,nbhaloanalyzed))
    If(.not.Allocated(haloCos_alpha))          Allocate(haloCos_alpha(nbhaloanalyzed))
    If(.not.Allocated(haloCenter_potential))   Allocate(haloCenter_potential(3,nbhaloanalyzed))
    If(.not.Allocated(haloCircular_values))    Allocate(haloCircular_values(4,nbhaloanalyzed))
    If(.not.Allocated(haloDelta_values))       Allocate(haloDelta_values(2,nbhaloanalyzed))
    If(.not.Allocated(haloRadial_profile))     Allocate(haloRadial_profile(NBIN1))
    If(.not.Allocated(haloRadial_profile178))  Allocate(haloRadial_profile178(NBIN2))
    If(.not.Allocated(haloVelocity_profile))   Allocate(haloVelocity_profile(NBIN1,nbhaloanalyzed))
    If(.not.Allocated(haloDensity_profile))    Allocate(haloDensity_profile(NBIN1,nbhaloanalyzed))
    If(.not.Allocated(haloDensity_profile178)) Allocate(haloDensity_profile178(NBIN2,nbhaloanalyzed))
    If(.not.Allocated(haloDensity_profile_integrated))    Allocate(haloDensity_profile_integrated(NBIN1,nbhaloanalyzed))
    If(.not.Allocated(haloDensity_profile_integrated178)) Allocate(haloDensity_profile_integrated178(NBIN2,nbhaloanalyzed))

    If(.not.Allocated(haloNpart))               Allocate(haloNpart(nbhaloanalyzed))
    If(.not.Allocated(haloVolume))              Allocate(haloVolume(nbhaloanalyzed))

    np = size(pos,2)

    if(do_periodic_boundary)then 
     !modify the halos cut by the periodic boundary conditions over x,y and z
     do j=1,3
       if((maxval(pos(j,:))-minval(pos(j,:)))>0.5) then
          do i=1,np
             if (pos(j,i)>0.5 ) then
                pos(j,i)=pos(j,i)-1.
             endif
          enddo
       endif
     enddo
    endif

    !Computation of halo's center-of-mass position and speed
    comx = 0
    comv = 0
    Do i = 1, np
       comx = comx + pos(:,i)
       comv = comv + vel(:,i)
    End Do

    comx = comx / real(np)
    comv = comv / real(np)



    if(do_periodic_boundary)then 
     !Shift to get mass_center between 0 and 1
     do j=1,3
        if(comx(j)<0.)then
           comx(j)  = comx(j)  + 1.
           pos(j,:) = pos(j,:) + 1.
        endif
     end do
    endif


    if(use_halo_finder_center_as_center_of_mass) then 
       !Will replace center
       !Find the currenthaloID in the hf_haloID catalog
       ihf=1
       do while(hf_haloID(ihf)/=currenthaloID.and.ihf<size(hf_halocompos,2))
          ihf=ihf+1
       end do

       if(hf_haloID(ihf)/=currenthaloID)then
          print*,'error cannot find currenthaloID in hf_haloID',ihf,hf_haloID(ihf),currenthaloID
          stop
       endif
       
       hfcomx=hf_halocompos(:,ihf)

       !Boundary condition
       if(do_periodic_boundary)then 
          do j=1,3
             if(hfcomx(j)-comx(j)>+0.5)hfcomx(j)=hfcomx(j)-1.
             if(hfcomx(j)-comx(j)<-0.5)hfcomx(j)=hfcomx(j)+1.
          end do
       endif
       

       !Warning: override comx
!       print*,'comx',comx
!      print*,'hfcomx',hfcomx
!       print*,' '
       comx=hfcomx


    endif



   

if(do_lightcone)then
!!!
!The following correction is needed if the halos belong to a lightcone,
!instead of a fixed redshift
    !Compute redshift of the current halo
    dtemp = sqrt(comx(1)**2 + comx(2)**2+ comx(3)**2) *box_len_Mpc !in Mpc/h
    zz   =  ZDIST(dtemp)

    !Compute the corrected units for the current halo
    simu_info%unit_l = simu_info%unit_l * (1./(1.+zz) /simu_info%aexp)
    simu_info%unit_d = simu_info%unit_d * (simu_info%aexp /(1./(1.+zz)))**3
    simu_info%unit_t = simu_info%unit_t * (1./(1.+zz) /simu_info%aexp)**2
    simu_info%aexp   = 1./(1.+zz)
    conv=(simu_info%h0/100./simu_info%aexp)
    Rho_m=simu_info%unit_d

    !Compute hubble expansion (in a^2 H/H0 units) at redshift z
    Hz = E_HUBBLE(simu_info%aexp) * simu_info%aexp**2
!!!
else !all halos belong to a fixed redshift (corresponding to the ramses snapshot)
    dtemp = -1.d0 !dtemp is not needed when halos belong to a fixed redshift
    zz    = 1./simu_info%aexp -1.
endif

    !Compute critical density (in g/cm^3) at redshift z
    Rho_critZ = (Hz / simu_info%aexp**2)**2 * Rho_crit0
    !print *,'z, Rho_m, Rho_critZ, Omega_m',zz, Rho_m, Rho_critZ, Rho_m/Rho_critZ

    M_vir= np*mass_part*simu_info%unit_d*simu_info%unit_l**3    !in grams
    R_vir=(3*M_vir/(4*PI*simu_info%delta_vir*Rho_m))**(1./3.)   !in cm
    V_vir = sqrt(G*M_vir/R_vir)                                 !in cm/s

    haloNpart(currenthalo)=np
    haloMvir(currenthalo) = real(M_vir/msun * simu_info%h0/100.,kind=8)      !in Msun/h
    haloRvir(currenthalo) = real(R_vir/(mpc/1d3)* simu_info%h0/100. ,kind=8)  !in kpc/h physical
    haloRvircom(currenthalo) = real(R_vir/(mpc/1d3)* simu_info%h0/100./simu_info%aexp ,kind=8)  !in kpc/h comoving
    haloVvir(currenthalo) = real(V_vir/1d5,kind=8)                           !in km/s

    haloIdentity(currenthalo) = currenthaloID
    haloCOMX(:,currenthalo) = real(comx*simu_info%unit_l/mpc *conv,kind=4) !in Mpc/h comoving
    haloCOMV(:,currenthalo) = real(comv*simu_info%unit_l/simu_info%unit_t/1d5,kind=4) !in km/s

    !Computation of R_max
    R_max=0
    R_tmp=0
    do i=1,np
       R_tmp=sqrt((pos(1,i)-comx(1))**2 +(pos(2,i)-comx(2))**2 +(pos(3,i)-comx(3))**2)
       if(R_tmp>R_max)R_max=R_tmp
    end do

    haloRmax(currenthalo) = real(R_max*simu_info%unit_l/mpc /(R_vir/mpc),kind=4) !physical/physical

    !Computation of DispX
    DispX=0
    do i=1,np
       Dispx = Dispx + (pos(1,i)-comx(1))**2 +(pos(2,i)-comx(2))**2 +(pos(3,i)-comx(3))**2
    end do

    DispX=sqrt(DispX/real(np)) !a better estimator consist in dividing by np-1
    haloDispX(currenthalo) = real(DispX*simu_info%unit_l/mpc /(R_vir/mpc),kind=4) !physical/physical

    !Computation of V_max
    V_max=0
    V_tmp=0
    do i=1,np
       V_tmp=sqrt((vel(1,i)-comv(1) + ((pos(1,i)-comx(1))*Hz))**2 & 
              & + (vel(2,i)-comv(2) + ((pos(2,i)-comx(2))*Hz))**2 &
              & + (vel(3,i)-comv(3) + ((pos(3,i)-comx(3))*Hz))**2)
       if(V_tmp>V_max)V_max=V_tmp
    end do

    haloVmax(currenthalo) = real(V_max*simu_info%unit_l/simu_info%unit_t/1d5/(V_vir/1d5),kind=4) ! km/s  / km/s

    !Computation of DispV
    DispV=0
    do i=1,np
       DispV = DispV + (vel(1,i)-comv(1) + ((pos(1,i)-comx(1))*Hz))**2 &
                   & + (vel(2,i)-comv(2) + ((pos(2,i)-comx(2))*Hz))**2 &
                   & + (vel(3,i)-comv(3) + ((pos(3,i)-comx(3))*Hz))**2
    end do

    DispV=sqrt(DispV/real(np)) !a better estimator consist in dividing by np-1
    haloDispV(currenthalo) = real(DispV*simu_info%unit_l/simu_info%unit_t/1d5 /(V_vir/1d5),kind=4) ! km/s  / km/s

    !Computation of total energy of halo
    Epot=0
    Ekin=0
    Etot=0
    Epot_max=0
    threshold=30000 !default 10 000
    threshold_calc=int(threshold/10.)

    Delta_X=(1./2.**simu_info%levelmin)*(1./2.**6.) !Smoothing length Big Assumption!!!

    !Computation of kinetic energy (including expansion)
    do i=1,np
       Ekin = Ekin + (vel(1,i)-comv(1)+(pos(1,i)-comx(1))*Hz)**2 &
                 & + (vel(2,i)-comv(2)+(pos(2,i)-comx(2))*Hz)**2 &
                 & + (vel(3,i)-comv(3)+(pos(3,i)-comx(3))*Hz)**2
    end do
    Ekin = 0.5*mass_part*Ekin            !1/2 * (mv²) ramses units
    Ekin = Ekin * simu_info%unit_d * (simu_info%unit_l**3) * (simu_info%unit_l**2) * (1./(simu_info%unit_t**2)) !cgs


    !Computation of the potential energy
    !Shankar modified the following, to find the particle with minimum potential in a halo
    if (np<threshold) then
       !if number of paticles in the halo is less, compute the potential energy using all particles
       do i=1,np
          Epot_i=0 !Epot_i is the potential energy of the i-th particle
          do j=1,np
            if (j /= i) then
               Epot_i = Epot_i + (mass_part**2)/sqrt((pos(1,i)-pos(1,j))**2 &
                    &   +(pos(2,i)-pos(2,j))**2 +(pos(3,i)-pos(3,j))**2 +Delta_X**2) !ramses units
            endif
          enddo

          Epot=Epot+Epot_i
          if (Epot_i > Epot_max) then !then this is the particle with minimum potential (so far) in this halo
             Epot_max = Epot_i
             center_pot(1:3)=pos(:,i) !store the location of the particle with minimum potential
          endif
       enddo
       Epot = -0.5*G*Epot * ((simu_info%unit_d * (simu_info%unit_l**3))**2) * (1./(simu_info%unit_l)) !G * (mm/|ri-rj|) m²/r is in cgs then mutiplied by G
    else
       !if number of particles in the halo is large, compute potential energy using a sample of particles
       if(.not.Allocated(rand_positions)) Allocate(rand_positions(3,threshold_calc),stat=ok)
       if(ok>0) then
          print*, procID,'Error in allocation of the rand_positions array'
          call MPI_FINALIZE(mpierr)
          stop
       endif
       if(.not.Allocated(random_index_arr)) Allocate(random_index_arr(threshold_calc),stat=ok)
       if(ok>0) then
          print*, procID,'Error in allocation of the random_index_arr array'
          call MPI_FINALIZE(mpierr)
          stop
       endif

       do i=1,threshold_calc
          test=.true.
          count_cpt=0
          if (count_cpt<10*threshold) then
             do while (test)
                call random_number(random)
                random_index=int(aint(np*random))+1
                test=.false.
                do j=1,i-1
                   if (random_index==random_index_arr(j)) then
                       test=.true.
                   endif
                enddo
                count_cpt=count_cpt+1
             enddo
          else
             test=.false.
             print*,procID,'The program could not find threshold_calc random particles through threshold.'
             call MPI_FINALIZE(mpierr)
             stop
          endif
          random_index_arr(i)=random_index
       enddo

       rand_positions(:,:)=pos(:,random_index_arr)

       do i=1,threshold_calc
          Epot_i=0 !Epot_i is the potential energy of the i-th particle
          do j=1,threshold_calc
            if (j /= i) then
               Epot_i=Epot_i+(mass_part**2)/sqrt((rand_positions(1,i)-rand_positions(1,j))**2 + &
                    & (rand_positions(2,i)-rand_positions(2,j))**2 + &
                    & (rand_positions(3,i)-rand_positions(3,j))**2 + Delta_X**2)
            endif
          enddo

          Epot=Epot+Epot_i
          if (Epot_i > Epot_max) then !then this is the particle with minimum potential (so far) in this halo
              Epot_max = Epot_i
              center_pot(1:3)=rand_positions(:,i) !store the location of the particle with minimum potential
          endif
       enddo
       deallocate(rand_positions,random_index_arr)

       !normalization to the right number of particles !!!CONVERT TO REAL*8 OTHERWISE BUG ABOVE 40000 PART RY 13/11/11
       Epot = - (real(np,kind=8)**2 - real(np,kind=8)) / (real(threshold_calc,kind=8)**2-real(threshold_calc,kind=8)) &
          &   * 0.5*G*Epot * ((simu_info%unit_d*(simu_info%unit_l**3))**2) * (1./(simu_info%unit_l))
    endif
    Etot=Ekin+Epot

    haloEpot(currenthalo) = real(Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)  !dimensionless
    haloEkin(currenthalo) = real(Ekin/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4) !dimensionless

    haloCenter_potential(:,currenthalo) = real(center_pot * simu_info%unit_l/mpc *conv,kind=4) !in Mpc/h comoving


   !Computation of halo density profile
   dndlnR=0
   dndlnR178=0
   density_profile178=0
   density_profile_integrated=0
   density_profile_integrated178=0
   V_circ=0
   V_circ_max=0
   M_circ_max=0
   R_circ_max=0
   Rho_circ_max=0
   Diff=1d10
   M178=0
   R178=0

   !set the central values of radial bins for profiles
   DO j=1,NBIN1
      haloRadial_profile(j) = real(exp(lmin_r+(j-0.5)*delta_lr)*(simu_info%h0/100.),kind=4) !in kpc/h physical 
   ENDDO

   DO j=1,NBIN2
      haloRadial_profile178(j) = real(exp(lmin_r178+(j-0.5)*delta_lr178),kind=4) !in R/R178
   ENDDO

   innermost_particles=0
   !Bin the halo particles in radial bins
   do i=1,np
      if(do_center_mass)then !halo density profile will be centered on center-of-mass
        R_tmp=sqrt((pos(1,i)-comx(1))**2+(pos(2,i)-comx(2))**2+(pos(3,i)-comx(3))**2)
      else  !halo density profile will be centered on minimum-of-potential
        R_tmp=sqrt((pos(1,i)-center_pot(1))**2+(pos(2,i)-center_pot(2))**2+(pos(3,i)-center_pot(3))**2)
      endif

      R_tmp = R_tmp *simu_info%unit_l /(mpc/1d3) !in kpc physical

     !First check if the particle is closer than even the smallest shell.
     if(log(R_tmp) .lt. lmin_r) then
        innermost_particles = innermost_particles +1
     else !Find the bin it belongs to.
        DO j=1,NBIN1
          if((lmin_r+(j-1)*delta_lr).le.log(R_tmp) .AND. log(R_tmp).lt.(lmin_r+j*delta_lr)) THEN
             dndlnR(j)=dndlnR(j)+1
          end if
        ENDDO
     endif
   end do

   !print*, currenthalo, np, innermost_particles, dndlnR, sum(dndlnR)

   !Calculate shell and integrated densities.
   DO j=1,NBIN1
      !Density in radial shells = Mass in shell / (4 pi R^2 dR)
      density_profile(j) = dndlnR(j) *mass_part*simu_info%unit_d*simu_info%unit_l**3 /msun *(simu_info%h0/100.) !Mass in radial shell
      !density_profile(j) = density_profile(j) /(4.*PI*(exp(lmin_r+(j-0.5)*delta_lr))**2) /delta_lr !Shell density [Msun/h / (kpc)^3]
      !instead of (4 pi R^2 dR), use the exact shell volume (4/3 pi (R2^3 - R1^3))
      density_profile(j) = density_profile(j) &
             & /(4./3.*PI*(exp(lmin_r+j*delta_lr)**3 - exp(lmin_r+(j-1)*delta_lr)**3)) !Shell density [Msun/h / (kpc)^3] !RY corrected bug

      density_profile(j) = density_profile(j) *msun /(simu_info%h0/100.) /(mpc/1d3)**3 !shell density in grams/cm3
      density_profile(j) = density_profile(j) /Rho_m     !shell density w.r.t background matter density at z
      !density_profile(j) = density_profile(j) /Rho_critZ !shell density w.r.t critical density at z

      !Spherically averaged density within radius R is M(<R)/(4/3 pi R^3)
      !Also calculate circular velocity profile, V_circ=sqrt(GM(<R)/R)
      !Also find the halo radius at which rho/rho_matter=178
      density_profile_integrated(j) = (innermost_particles+sum(dndlnR(1:j))) *mass_part*simu_info%unit_d*simu_info%unit_l**3 !M(<R) in grams
      R_tmp = exp(lmin_r+j*delta_lr) *(mpc/1d3) !in cm physical

      V_circ(j) = sqrt(G *density_profile_integrated(j) /R_tmp) !V_circ in cm/sec
      if(V_circ(j) > V_circ_max)then
        V_circ_max = V_circ(j) !V_circ_max in cm/s
        M_circ_max = density_profile_integrated(j) /msun *(simu_info%h0/100.) !M(<R) in Msun/h
        R_circ_max = R_tmp /(mpc/1d3) *(simu_info%h0/100.) !in kpc/h physical
        Rho_circ_max = density_profile_integrated(j) /(4./3.*PI*R_tmp**3) !in grams/cm3
        Rho_circ_max = Rho_circ_max /Rho_m !w.r.t background matter density at z
      endif

      density_profile_integrated(j) = density_profile_integrated(j) /(4./3.*PI*R_tmp**3) !in grams/cm3
      density_profile_integrated(j) = density_profile_integrated(j) /Rho_m !w.r.t background matter density at z
      !density_profile_integrated(j) = density_profile_integrated(j) /Rho_critZ !w.r.t critical density at z
      if(abs(density_profile_integrated(j)-Delta_matter) < Diff) THEN
         R178 = R_tmp /(mpc/1d3) !in kpc physical
         M178 = (4./3.*PI*R_tmp**3) *density_profile_integrated(j) *Rho_m  !M(<R178) in grams
         M178 = M178  /msun *(simu_info%h0/100.) !M(<R178) in Msun/h
         Diff=abs(density_profile_integrated(j)-Delta_matter)
      endif
   ENDDO
   
   !Interpolate R178
   jstart=10
   j=1
   neverabovedelta=.true.
   do while ((density_profile_integrated(j)>Delta_matter.and.j<NBIN1).or.(neverabovedelta.and.j<NBIN1).or.(j<jstart.and.j<NBIN1) )
      if(neverabovedelta)then
         if(density_profile_integrated(j)>Delta_matter) then
            neverabovedelta=.false.
         endif
      endif
      j=j+1
   end do
   if (j==NBIN1) then 
      print*,'WARNING Rdelta not found in range of radius up to bin',NBIN1, density_profile_integrated(j),Delta_matter,currenthaloID
      print*,'Will extrapolate'
   endif
   if(j==jstart) then
      print*,'WARNING Rdelta already found at first bin',jstart, density_profile_integrated(j),Delta_matter,currenthaloID
      print*,'Will extrapolate'   
   endif
   
   if(Delta_matter<0..or.density_profile_integrated(j-1)<0..or.density_profile_integrated(j)<0.)then
      print*,'erreur sign density',Delta_matter,density_profile_integrated(j-1),density_profile_integrated(j),currenthaloID
      stop
   endif

   R_tmp=(lmin_r+j*delta_lr)*    (dlog(Delta_matter)-dlog(density_profile_integrated(j-1)))/(dlog(density_profile_integrated(j))  -dlog(density_profile_integrated(j-1)))+ &
        &(lmin_r+(j-1)*delta_lr)*(dlog(Delta_matter)-dlog(density_profile_integrated(j)))   /(dlog(density_profile_integrated(j-1))-dlog(density_profile_integrated(j)))
   R_tmp=exp(R_tmp)*(mpc/1d3) !in cm physical
   R178= R_tmp /(mpc/1d3) !in kpc physical
   M178 = (4./3.*PI*R_tmp**3) *Delta_matter *Rho_m 
   M178 = M178  /msun *(simu_info%h0/100.)


   !write(*,*)'M178,R178',currenthalo,currenthaloID,M178,R178

   haloVelocity_profile(:,currenthalo) = real(V_circ/1d5,kind=4)     !in km/s

   haloCircular_values(1,currenthalo) = real(M_circ_max,kind=4)      !in Msun/h
   haloCircular_values(2,currenthalo) = real(R_circ_max,kind=4)      !in kpc/h physical
   haloCircular_values(3,currenthalo) = real(V_circ_max/1d5,kind=4)  !in km/s
   haloCircular_values(4,currenthalo) = real(Rho_circ_max,kind=4)    !in background matter

   haloDelta_values(1,currenthalo) = real(M178,kind=4) !in Msun/h
   haloDelta_values(2,currenthalo) = real(R178*simu_info%h0/100.,kind=4) !in kpc/h physical

   innermost_particles=0
   !Now, re-bin the halo particles in radial bins. This time, bins are in r/r178 units.
   do i=1,np
      if(do_center_mass)then !halo density profile will be centered on center-of-mass
        R_tmp=sqrt((pos(1,i)-comx(1))**2+(pos(2,i)-comx(2))**2+(pos(3,i)-comx(3))**2)
      else  !halo density profile will be centered on minimum-of-potential
        R_tmp=sqrt((pos(1,i)-center_pot(1))**2+(pos(2,i)-center_pot(2))**2+(pos(3,i)-center_pot(3))**2)
      endif

      R_tmp = R_tmp *simu_info%unit_l /(mpc/1d3) !in kpc physical
      R_tmp = R_tmp /R178 !in R/R178 units

     !First check if the particle is closer than even the smallest shell.
     if(log(R_tmp) .lt. lmin_r178) then
        innermost_particles = innermost_particles +1
     else !Find the bin it belongs to.
        DO j=1,NBIN2
          if((lmin_r178+(j-1)*delta_lr178).le.log(R_tmp) .AND. log(R_tmp).lt.(lmin_r178+j*delta_lr178)) THEN
             dndlnR178(j)=dndlnR178(j)+1
          end if
        ENDDO
     endif
   end do

   !print*, currenthalo, np, innermost_particles, dndlnR178, sum(dndlnR178)

   !Calculate shell and integrated densities. Radial bins are in R/R178 units.
   DO j=1,NBIN2
      !Density in radial shells = Mass in shell / (4 pi R^2 dR)
      density_profile178(j) = dndlnR178(j) *mass_part*simu_info%unit_d*simu_info%unit_l**3 /msun *(simu_info%h0/100.) !Mass in radial shell [Msun/h]
      !density_profile178(j) = density_profile178(j) /(4.*PI*(exp(lmin_r178+(j-0.5)*delta_lr178))**3) /delta_lr178 /R178**3 !Shell density [Msun/h / (kpc)^3]
      !instead of (4 pi R^2 dR), use the exact shell volume (4/3 pi (R2^3 - R1^3))
      density_profile178(j) = density_profile178(j) /R178**3 &
             & /(4./3.*PI*(exp(lmin_r178+j*delta_lr178)**3 - exp(lmin_r178+(j-1)*delta_lr178)**3))  !Shell density [Msun/h / (kpc)^3]

      density_profile178(j) = density_profile178(j) *msun /(simu_info%h0/100.) /(mpc/1d3)**3 !shell density in grams/cm3
      density_profile178(j) = density_profile178(j) /Rho_m     !shell density w.r.t background matter density at z
      !density_profile178(j) = density_profile178(j) /Rho_critZ !shell density w.r.t critical density at z

      !Spherically averaged density within radius R is M(<R)/(4/3 pi R^3)
      density_profile_integrated178(j) = (innermost_particles+sum(dndlnR178(1:j))) *mass_part*simu_info%unit_d*simu_info%unit_l**3 !M(<R) in grams
      R_tmp = exp(lmin_r178+j*delta_lr178) *R178 *(mpc/1d3) !physical units (in centimeters)

      density_profile_integrated178(j) = density_profile_integrated178(j) /(4./3.*PI*R_tmp**3) !in grams/cm3
      density_profile_integrated178(j) = density_profile_integrated178(j) /Rho_m     !w.r.t background matter density at z
      !density_profile_integrated178(j) = density_profile_integrated178(j) /Rho_critZ !w.r.t critical density at z
   ENDDO

   !Computation of halo volume
   compute_volume=.false. !Set to false if you want to avoid computing volume which is memory and cpu intensive

   if(compute_volume)then
      !Rough estimate of the volume
      bfof=0.2*(200./Delta_matter)**3 !bfof is the size of the grid to compute volume. Here, assume a very naive law to pave optimally the volume.
      bfof=bfof*mass_part**(1./3.)
      volume=cmp_volume_fofgroup(np,pos,bfof,.false.) !everything should be in ramses units (ie box unit)
      haloVolume(currenthalo) = real(volume*(simu_info%unit_l/mpc)**3/(4./3.*PI*(R_vir/mpc)**3) , kind=8) !normalized by 4/3*pi*normalisation_radius^3
   else
      haloVolume(currenthalo)=0.
   endif

   haloDensity_profile(:,currenthalo)               = real(density_profile,kind=4)
   haloDensity_profile178(:,currenthalo)            = real(density_profile178,kind=4)
   haloDensity_profile_integrated(:,currenthalo)    = real(density_profile_integrated,kind=4)
   haloDensity_profile_integrated178(:,currenthalo) = real(density_profile_integrated178,kind=4)


   


   !Computation of halo's angular momentum
   ang_momentum=0
   do i=1,np
      !Computation of moment over x, y and z
      ang_momentum(1) = ang_momentum(1) + &
           & mass_part*(pos(2,i) - comx(2)) * (vel(3,i) - comv(3)) - &
           & mass_part*(pos(3,i) - comx(3)) * (vel(2,i) - comv(2))

      ang_momentum(2) = ang_momentum(2) + &
           & mass_part*(pos(3,i) - comx(3)) * (vel(1,i) - comv(1)) - &
           & mass_part*(pos(1,i) - comx(1)) * (vel(3,i) - comv(3))

      ang_momentum(3) = ang_momentum(3) + &
           & mass_part*(pos(1,i) - comx(1)) * (vel(2,i) - comv(2)) - &
           & mass_part*(pos(2,i) - comx(2)) * (vel(1,i) - comv(1))
   end do

   !physical renormalisation in cgs of the angular momentum
   ang_momentum = ang_momentum * simu_info%unit_d * (simu_info%unit_l**5)/simu_info%unit_t !cgs

   jtot=sqrt(ang_momentum(1)**2 + ang_momentum(2)**2 + ang_momentum(3)**2)

   !normalisation of J in lambda : Lamda=J*sqrt(E)/(G*M_vir**(5/2))
   Lambda = ang_momentum * sqrt(abs(Etot))/(G*M_vir**(5./2.))
   haloLambda(:,currenthalo) = real(Lambda,kind=4)

   !normalisation of J in lambda' : Lamda'=J/(sqrt(2)*M_vir*V_vir*R_vir)
   Lambda_prime = ang_momentum / (sqrt(2.)*M_vir*R_vir*V_vir)
   haloLambda_prime(:,currenthalo) = real(Lambda_prime,kind=4)

   Inert_Tens(:,:)=0
   ! Computation of halo's inertial tensor
   do i=1,np
      Inert_Tens(1,1)=Inert_Tens(1,1)+ (pos(1,i)-comx(1)) * (pos(1,i)-comx(1))
      Inert_Tens(1,2)=Inert_Tens(1,2)+ (pos(1,i)-comx(1)) * (pos(2,i)-comx(2))
      Inert_Tens(1,3)=Inert_Tens(1,3)+ (pos(1,i)-comx(1)) * (pos(3,i)-comx(3))
      Inert_Tens(2,1)=Inert_Tens(2,1)+ (pos(2,i)-comx(2)) * (pos(1,i)-comx(1))
      Inert_Tens(2,2)=Inert_Tens(2,2)+ (pos(2,i)-comx(2)) * (pos(2,i)-comx(2))
      Inert_Tens(2,3)=Inert_Tens(2,3)+ (pos(2,i)-comx(2)) * (pos(3,i)-comx(3))
      Inert_Tens(3,1)=Inert_Tens(3,1)+ (pos(3,i)-comx(3)) * (pos(1,i)-comx(1))
      Inert_Tens(3,2)=Inert_Tens(3,2)+ (pos(3,i)-comx(3)) * (pos(2,i)-comx(2))
      Inert_Tens(3,3)=Inert_Tens(3,3)+ (pos(3,i)-comx(3)) * (pos(3,i)-comx(3))
   end do
   inert_tmp=Inert_Tens/real(np)
   !search the eigen values and vectors using the subroutines diag and eigsrt
   call diag(Inert_tmp,3,3,Ein_Val_Inert_Tens,Ein_Vec_Inert_Tens,nrot) !!!Yann divided by mass_halo!!!
   call eigsrt(Ein_Val_Inert_Tens,Ein_Vec_Inert_Tens,3,3)

   Ein_Val_Inert_Tens=sqrt(Ein_Val_Inert_Tens)

   haloEin_Val_Inert_Tens(:,currenthalo) = real(Ein_Val_Inert_Tens *simu_info%unit_l/mpc  /(R_vir/mpc),kind=4)
   haloEin_Vec_Inert_Tens(:,currenthalo) = real(Ein_Vec_Inert_Tens(:,3),kind=4)


   cos_alpha=(Ein_Vec_Inert_Tens(1,3)*ang_momentum(1) + &
        & Ein_Vec_Inert_Tens(2,3)*ang_momentum(2) + &
        & Ein_Vec_Inert_Tens(3,3)*ang_momentum(3) ) &
        & /sqrt(( Ein_Vec_Inert_Tens(1,3)**2 + Ein_Vec_Inert_Tens(2,3)**2 + Ein_Vec_Inert_Tens(3,3)**2 ) * &
        & ( ang_momentum(1)**2 + ang_momentum(2)**2 + ang_momentum(3)**2 ))

   haloCos_alpha(currenthalo) = cos_alpha

!!!Write single ascii file with everything
!! ID parallel FoF, number of particles
!! distance (Mpc/h comoving), redshift
!! Mvir (Msun/h), Rvir (kpc/h physical), Vvir (km/s) (Rvir and Vvir from Mvir assuming Deltavir)
!! center-of-mass x y z (Mpc/h comoving)
!! vx vy vz (km/s)
!! R_max/R_vir V_max/V_vir !physical/physical
!! sigmar/Rvir sigmav/V_vir !physical/physical
!! Ekin/(1/2 Mvir Vvir^2)   -Epot/(1/2 Mvir Vvir^2) !physical/physical
!! lambda_prime, lambda_primex,lambda_primey, lambda_primez !physical/physical
!! a/Rvir, b/Rvir, c/Rvir !physical/physical
!! minorx,minory,minorz
!! minimum-of-potential x y z (Mpc/h comoving)
!! M_circ_max (Msun/h), R_circ_max (kpc/h physical), V_circ_max (km/s), Rho_circ_max (background)
!! M178 (Msun/h), R178 (kpc/h physical)

   if(ascii) then
      open(400,file='halo_properties.dat',form='formatted',status='unknown',access='append')
      open(401,file='haloDensity_profile.dat',form='formatted',status='unknown',access='append')
      open(402,file='haloDensity_profile_integrated.dat',form='formatted',status='unknown',access='append')
      open(403,file='haloDensity_profile178.dat',form='formatted',status='unknown',access='append')
      open(404,file='haloDensity_profile_integrated178.dat',form='formatted',status='unknown',access='append')
      open(405,file='haloVelocity_profile.dat',form='formatted',status='unknown',access='append')
      
      if(procID==0 .AND. currenthalo==1)then !write header
         write(400,'(38a18)')'idp','npart','Dist','z','mass','r178','v178','x','y','z','vx','vy','vz','rmax','vmax'&
              &,'sigmapos','sigmavel','ekin','epot','jtot','jx','jy','jz','a','b','c'&
              &,'minorx','minory','minorz','potential_x','potential_y','potential_z'&
              &,'M_circ_max','R_circ_max','V_circ_max','Rho_circ_max','M178','R178'
      endif
      
      write(400,'((2i15,10000g18.7))') currenthaloID, np&
           &,real(dtemp,kind=4), real(zz,kind=4)&
           &,real(M_vir/msun *(simu_info%h0/100.),kind=4)&
           &,real(R_vir/(mpc/1d3) *(simu_info%h0/100.),kind=4)&
           &,real(V_vir/1d5,kind=4)&
           &,real(comx *simu_info%unit_l/mpc *conv,kind=4)&
           &,real(comv *simu_info%unit_l/simu_info%unit_t/1d5,kind=4)&
           &,real(R_max*simu_info%unit_l/mpc /(R_vir/mpc),kind=4)&
           &,real(V_max*simu_info%unit_l/simu_info%unit_t/1d5 /(V_vir/1d5),kind=4)&
           &,real(DispX*simu_info%unit_l/mpc /(R_vir/mpc),kind=4)&
           &,real(DispV*simu_info%unit_l/simu_info%unit_t/1d5 /(V_vir/1d5),kind=4)&
           &,real(Ekin/(msun*1d10)  /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(-Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(jtot/(sqrt(2.)*M_vir*R_vir*V_vir),kind=4)&
           &,real(ang_momentum /(sqrt(2.)*M_vir*R_vir*V_vir),kind=4)&
           &,real(Ein_Val_Inert_Tens *simu_info%unit_l/mpc /(R_vir/mpc),kind=4)&
           &,real(Ein_Vec_Inert_Tens(:,3),kind=4)&
           &,real(center_pot *simu_info%unit_l/mpc *conv,kind=4)&
           &,real(M_circ_max,kind=4), real(R_circ_max,kind=4)&
           &,real(V_circ_max/1d5,kind=4), real(Rho_circ_max,kind=4)&
           &,real(M178,kind=4), real(R178 *(simu_info%h0/100.),kind=4)
      close(400)
      
      write(401,'((2i15,10000g18.7))') currenthaloID, np&
           &,real(dtemp,kind=4), real(zz,kind=4)&
           &,real(M_vir/msun *(simu_info%h0/100.),kind=4)&
           &,real(M178,kind=4), real(R178 *(simu_info%h0/100.),kind=4)&
           &,real(Ekin/(msun*1d10)  /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(-Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(jtot/(sqrt(2.)*M_vir*R_vir*V_vir),kind=4), real(density_profile,kind=4)
      close(401)
      
      write(402,'((2i15,10000g18.7))') currenthaloID, np&
           &,real(dtemp,kind=4), real(zz,kind=4)&
           &,real(M_vir/msun *(simu_info%h0/100.),kind=4)&
           &,real(M178,kind=4), real(R178 *(simu_info%h0/100.),kind=4)&
           &,real(Ekin/(msun*1d10)  /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(-Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(jtot/(sqrt(2.)*M_vir*R_vir*V_vir),kind=4), real(density_profile_integrated,kind=4)
      close(402)
      
      write(403,'((2i15,10000g18.7))') currenthaloID, np&
           &,real(dtemp,kind=4), real(zz,kind=4)&
           &,real(M_vir/msun *(simu_info%h0/100.),kind=4)&
           &,real(M178,kind=4), real(R178 *(simu_info%h0/100.),kind=4)&
           &,real(Ekin/(msun*1d10)  /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(-Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(jtot/(sqrt(2.)*M_vir*R_vir*V_vir),kind=4), real(density_profile178,kind=4)
      close(403)
      
      write(404,'((2i15,10000g18.7))') currenthaloID, np&
           &,real(dtemp,kind=4), real(zz,kind=4)&
           &,real(M_vir/msun *(simu_info%h0/100.),kind=4)&
           &,real(M178,kind=4), real(R178 *(simu_info%h0/100.),kind=4)&
           &,real(Ekin/(msun*1d10)  /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(-Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(jtot/(sqrt(2.)*M_vir*R_vir*V_vir),kind=4), real(density_profile_integrated178,kind=4)
      close(404)
      
      write(405,'((2i15,10000g18.7))') currenthaloID, np&
           &,real(dtemp,kind=4), real(zz,kind=4)&
           &,real(M_vir/msun *(simu_info%h0/100.),kind=4)&
           &,real(M178,kind=4), real(R178 *(simu_info%h0/100.),kind=4)&
           &,real(Ekin/(msun*1d10)  /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(-Epot/(msun*1d10) /(1./2.*M_vir/msun*(V_vir/1d5)**2),kind=4)&
           &,real(jtot/(sqrt(2.)*M_vir*R_vir*V_vir),kind=4), real(V_circ/1d5,kind=4)
      close(405)
   endif
   
  End Subroutine halo_properties
  !=======================================================================


  subroutine eigsrt(d,v,n,np)
    integer :: n,np
    real (kind=8) :: p, d(np),v(np,np)
    integer :: i,j,k

    do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
            if(d(j) .ge. p) then
                k=j
                p=d(j)
            endif
        enddo
        if (k .ne. i) then
            d(k)=d(i)
            d(i)=p
            do j=1,n
                p=v(j,i)
                v(j,i)=v(j,k)
                v(j,k)=p
            enddo
        endif
    enddo
    return
  end subroutine eigsrt
  !=======================================================================

  subroutine diag(a,n,np,d,v,nrot)
    !compute the eigen values 'd' and eigen vectors 'v' of 'a'
    integer :: n,np,nrot,NMAX
    real (kind=8) ::  a(np,np),d(np),v(np,np)
    parameter (NMAX=500)
    integer :: i,ip,iq,j
    real (kind=8) :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

    do ip=1,n
        do iq=1,n
            v(ip,iq)=0.
        enddo
        v(ip,ip)=1.
    enddo
    do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.

    enddo
    nrot=0
    do i=1,50
        sm=0.
        do ip=1,n-1
            do iq=ip+1,n
                sm=sm+abs(a(ip,iq))
            enddo
        enddo
        if(sm.eq.0.)return
            if(i.lt.4)then
                tresh=0.2*sm/n**2
            else
                tresh=0.
            endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                    if(abs(h)+g.eq.abs(h))then
                        t=a(ip,iq)/h
                    else
                        theta=0.5*h/a(ip,iq)
                        t=1./(abs(theta)+sqrt(1.+theta**2))
                        if(theta.lt.0.)t=-t
                    endif
                    c=1./sqrt(1+t**2)
                    s=t*c
                    tau=s/(1.+c)
                    h=t*a(ip,iq)
                    z(ip)=z(ip)-h
                    z(iq)=z(iq)+h
                    d(ip)=d(ip)-h
                    d(iq)=d(iq)+h
                    a(ip,iq)=0.
                    do j=1,ip-1
                        g=a(j,ip)
                        h=a(j,iq)
                        a(j,ip)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
                    enddo
                    do j=ip+1,iq-1
                        g=a(ip,j)
                        h=a(j,iq)
                        a(ip,j)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
                    enddo
                    do j=iq+1,n
                        g=a(ip,j)
                        h=a(iq,j)
                        a(ip,j)=g-s*(h+g*tau)
                        a(iq,j)=h+s*(g-h*tau)
                    enddo
                    do j=1,n
                        g=v(j,ip)
                        h=v(j,iq)
                        v(j,ip)=g-s*(h+g*tau)
                        v(j,iq)=h+s*(g-h*tau)
                    enddo
                    nrot=nrot+1
                endif
            enddo
        enddo
        do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.
        enddo
    enddo
   !pause 'too many iterations in jacobi'
   return
  end subroutine diag
  !=======================================================================


  real(kind=8) function cmp_volume_fofgroup(npart,pos,bfof,verbose)
  !Compute the volume encompassed by a distribution of npart particles at
  !position pos.
  !The volume will be computed by considering a grid of size dx=dy=dz=bfof where
  !bfof is the typical mean 
  !interparticle seperation. verbose activates verbose mode.
  !The method is from More et al, 2011 (with a very small variation): consider
  !the particles projected on a 
  !2D mesh (let's say in x-y plane). Then for each 2D cells compute the
  !extension of the particles
  !along the third dimension (zmax-zmin). The volume is the sum of the volume of
  !all 1D pencils dx*dy*(dz*(izmax-izmin)).
  !The projection is done on the 3 axes and the maximum of the 3 volume is
  !taken. This methods minimize
  !underestimates (due to Poisson fluctuation) when one computes the sum of the
  !volume of all non-empty cells.
  !There are still some fluctuations on the surface of the volume: better
  !solution could be for instance concave hull computation, or 
  !tessellation or some monte-carlo methods maybe?
  !Y.Rasera 07/09/2015

  implicit none
  !input
  integer::npart   !nb of particles
  real,dimension(3,npart)::pos  !input arrays first dimension are the dimensions second are part number
  real(kind=8)::bfof  !linking length which will set the typical grid size
  logical::verbose   !verbose mode

  !internal variables
  real::xminini,xmaxini,yminini,ymaxini,zminini,zmaxini  !min and max positions
  real::xmin,xmax,ymin,ymax,zmin,zmax !min and max effectively used
  integer,dimension(:,:,:),allocatable::cube !0 if cell is empty 1 if cell is non-empty
  integer::idim,ix,iy,iz,nx,ny,nz,imin,imax !index in each dimensions
  integer:: i !index of particles
  real::dx,dy,dz !size of cells
  real::ratiox,ratioy,ratioz ! to 
  real::restx,resty,restz !free space between edge of grid and last particles:
  !to center the grid so as to encompass the particles distribution and let
  !equal space on each side restxyz/2
  real::vol,vol0,vol1,vol2,vol3 ! volume estimates
  logical::mem_security
  integer::size_security
  real::fact
  logical::saturation,continueafter

  if(verbose) then
     print*,'npart=',npart
     print*,'bfof =',bfof
  endif

  cmp_volume_fofgroup=0.

  !Compute minimum/maximum of positions
  xminini=minval(pos(1,:))
  xmaxini=maxval(pos(1,:))
  yminini=minval(pos(2,:))
  ymaxini=maxval(pos(2,:))
  zminini=minval(pos(3,:))
  zmaxini=maxval(pos(3,:))

  if(verbose)then
     print*,'min(posx) max(posx)=',xminini,xmaxini
     print*,'min(posy) max(posy)=',yminini,ymaxini
     print*,'min(posz) max(posz)=',zminini,zmaxini
  endif

  !Set grid spacing to bfof
  dx=bfof
  dy=bfof
  dz=bfof
  if(verbose)then
     print*,'dx=',dx
     print*,'dy=',dy
     print*,'dz=',dz
  endif

  !Compute number of grids using ceiling so that grids are beyond particles.
  ratiox=(xmaxini-xminini)/bfof
  ratioy=(ymaxini-yminini)/bfof
  ratioz=(zmaxini-zminini)/bfof
  nx=ceiling(ratiox)
  ny=ceiling(ratioy)
  nz=ceiling(ratioz)
  if(verbose)then
     print*,'nx=',nx
     print*,'ny=',ny
     print*,'nz=',nz
  endif
  
  !Security
  continueafter=.true.
  mem_security=.true.
  saturation=.false.
  size_security=100
  if(mem_security) then
     if(nx*ny*nz.gt.size_security*npart)then
        print*,'WARNING SIZE ARRAY TO COMPUTE VOLUME IS TOO LARGE',nx*ny*nz,size_security*npart,npart,size_security
        if(saturation) then
           print*,'WILL SATURATE TO ROUGHLY',size_security*npart 
           
           fact=real(nx*ny*nz)/real(size_security*npart)
           
           bfof=bfof/fact**(1./3.)
           if(verbose) then
              print*,'NEW bfof IS =',bfof
           endif
           
           !Set grid spacing to bfof
           dx=bfof
           dy=bfof
           dz=bfof
           if(verbose)then
              print*,'dx=',dx
              print*,'dy=',dy
              print*,'dz=',dz
           endif
           
           !Compute number of grids using ceiling so that grids are beyond particles.
           ratiox=(xmaxini-xminini)/bfof
           ratioy=(ymaxini-yminini)/bfof
           ratioz=(zmaxini-zminini)/bfof
           nx=ceiling(ratiox)
           ny=ceiling(ratioy)
           nz=ceiling(ratioz)
           if(verbose)then
              print*,'nx=',nx
              print*,'ny=',ny
              print*,'nz=',nz
           endif
        else
           print*,'will not compute volume for this particular halo'
           continueafter=.false.
        endif
     endif
  endif

  if(continueafter) then
     !Allocate cube array (recentring requires +1)
     allocate(cube(nx+1,ny+1,nz+1))
     cube=0
     
     !Compute the empty space between edge of last grid and last particles
     !Then recentered to let equal space on each side of the particles distribution
     restx=nx*bfof-(xmaxini-xminini)
     resty=ny*bfof-(ymaxini-yminini)
     restz=nz*bfof-(zmaxini-zminini)
     if(verbose)then
        print*,'restx',restx
        print*,'resty',resty
        print*,'restz',restz
     endif
     xmin=xminini-restx/2.
     xmax=xmaxini+restx/2.
     ymin=yminini-resty/2.
     ymax=ymaxini+resty/2.
     zmin=zminini-restz/2.
     zmax=zmaxini+restz/2.
     if(verbose)then
        print*,'xminxmax',xmin,xmax
        print*,'yminymax',ymin,ymax
        print*,'zminzmax',zmin,zmax
     endif
     
     !Set cube to 1 in each non-empty cells
     do i=1,npart
        ix=int((pos(1,i)-xmin)/dx)+1
        iy=int((pos(2,i)-ymin)/dy)+1
        iz=int((pos(3,i)-zmin)/dz)+1     
        cube(ix,iy,iz)=1
     end do
     
     !Naive estimate of the volume by counting non-empty cells
     vol0=sum(cube)*dx*dy*dz
     if(verbose)then
        print*,'naive volume estimate',vol0
     endif

     !Compute the volume by projection along z, computation of extension along z,
     !and sum of dx*dy*extensionz
     vol1=0.
     do ix=1,nx+1
        do iy=1,ny+1
           
           imin=nz+2
           do iz=1,nz+1
              idim=nz-iz+2
              if((cube(ix,iy,idim)==1))imin=idim
           end do
           imax=0
           do iz=1,nz+1
              idim=iz
              if((cube(ix,iy,idim)==1))imax=idim
           end do
           if(imin<=imax)vol1=vol1+(imax-imin+1)
           
        end do
     end do
     
     !Compute the volume by projection along y, computation of extension along y,
     !and sum of dx*dz*extensiony
     vol2=0.
     do ix=1,nx+1
        do iz=1,nz+1
           
           imin=ny+2
           do iy=1,ny+1
              idim=ny-iy+2
              if((cube(ix,idim,iz)==1))imin=idim
           end do
           imax=0
           do iy=1,ny+1
              idim=iy
              if((cube(ix,idim,iz)==1))imax=idim
           end do
           if(imin<=imax)vol2=vol2+(imax-imin+1)
           
        end do
     end do
     
     !Compute the volume by projection along x, computation of extension along x,
     !and sum of dy*dz*extensionx
     vol3=0.
     do iy=1,ny+1
        do iz=1,nz+1
           
           imin=nx+2
           do ix=1,nx+1
              idim=nx-ix+2
              if((cube(idim,iy,iz)==1))imin=idim
           end do
           imax=0
           do ix=1,nx+1
              idim=ix
              if((cube(idim,iy,iz)==1))imax=idim
           end do
           if(imin<=imax)vol3=vol3+(imax-imin+1)
           
        end do
     end do
     
     !Take the max of the 3 volumes (in grid units)
     vol=maxval((/vol1,vol2,vol3/))
     
     !Final volume (in input units)
     cmp_volume_fofgroup=vol*dx*dy*dz
     
     !Deallocate cube
     deallocate(cube)

     if(verbose)then
        print*,'volume estimates',vol1*dx*dy*dz,vol2*dx*dy*dz,vol3*dx*dy*dz
     endif
  endif

end function cmp_volume_fofgroup




End Module modfunctions


