pro test
rd_pfof_halo_snaporcone_part_prop_hdf5,ppf,pff,file='/obs/yrasera/proj/fulluniverse/code/haloanalyzer_new/pfof_halo_snap_part_prop_boxlen82.03125_n128_lcdmw7.h5'
rd_pfof_halo_snaporcone_part_prop_hdf5,pps,pfs,file='/obs/yrasera/proj/fulluniverse/code/haloanalyzer_new/psod_halo_snap_part_prop_boxlen82.03125_n128_lcdmw7.h5'

rd_pfof_halo_snaporcone_part_data_hdf5,pf,file='/data/yrasera/prepa4096/boxlen82.03125_n128_lcdmw7/post/fof/output_00011/pfof_halo_snap_part_data_boxlen82.03125_n128_lcdmw7_00000.h5',icpu=0L,ncpu=8L,/id,/pos,/vel
rd_pfof_halo_snaporcone_part_data_hdf5,ps,file='/data/yrasera/prepa4096/boxlen82.03125_n128_lcdmw7/post/sod/test_00011/psod_halo_snap_part_data_boxlen82.03125_n128_lcdmw7_00000.h5',icpu=0L,ncpu=8L,/id,/pos,/vel

;rd_pfof_halo_snaporcone_part_hfprop_hdf5,hfps,file='/data/yrasera/prepa4096/boxlen82.03125_n128_lcdmw7/post/sod/test_00011/psod_halo_snap_part_hfprop_boxlen82.03125_n128_lcdmw7.h5

;SOD VS FOF
plot,pps.POSITION_CENTER_OF_MASS_HALO(0,*),pps.POSITION_CENTER_OF_MASS_HALO(1,*),psym=3,xr=[0,82.03125],yr=[0,82.03125],/xs,/ys
oplot,ppf.POSITION_CENTER_OF_MASS_HALO(0,*),ppf.POSITION_CENTER_OF_MASS_HALO(1,*),color=2,psym=3

his,alog10(pps.mass_halo),/ylog
his,alog10(ppf.mass_halo),/ylog,/opl,color=2

stop

numhalof=0L
ihalof=ppf.identity_halo(numhalof)
indpf=where(pf.idfofhalo eq ihalof(0),nokf)
plot,pf.xp(indpf,0),pf.xp(indpf,1),psym=3,/xs,/ys
print,'minmax(x)',minmax(pf.xp(indpf,0))
print,'minmax(y)',minmax(pf.xp(indpf,1))
print,'minmax(z)',minmax(pf.xp(indpf,2))

stop
dfs2=(pps.POSITION_CENTER_OF_MASS_HALO(0,*)-ppf.POSITION_CENTER_OF_MASS_HALO(0,numhalof))^2+(pps.POSITION_CENTER_OF_MASS_HALO(1,*)-ppf.POSITION_CENTER_OF_MASS_HALO(1,numhalof))^2+(pps.POSITION_CENTER_OF_MASS_HALO(2,*)-ppf.POSITION_CENTER_OF_MASS_HALO(2,numhalof))^2
numhalos=where(dfs2 eq min(dfs2))
ihalos=pps.identity_halo(numhalos)
indps=where(ps.idfofhalo eq ihalos(0))
oplot,ps.xp(indps,0),ps.xp(indps,1),psym=3,color=2

stop

;CHECK PROP

print,ppf.mass_halo(numhalof)/min(ppf.mass_halo)*100,nokf,ppf.number_particles_halo(numhalof)
print,ppf.NORMALISATION_MASS_HALO(numhalof)/min(ppf.NORMALISATION_MASS_HALO)*100,nokf
print,ppf.NORMALISATION_RADIUS_PHYSICAL_HALO(numhalof),(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.)*3.17682715523738E26/3.08568E21*0.72
print,ppf.NORMALISATION_RADIUS_COMOVING_HALO(numhalof),(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.)*82.03125d0*1000.

print,ppf.NORMALISATION_VELOCITY_HALO(numhalof),sqrt(6.673099761655976E-8*(nokf/128.^3*3.39840378010395d-30*3.17682715523738d26^3)/(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.)/3.17682715523738E26)/1e5
print,ppf.POSITION_CENTER_OF_MASS_HALO(0,numhalof),mean(pf.xp(indpf,0),/double)*82.03125d0
print,ppf.POSITION_CENTER_OF_MASS_HALO(1,numhalof),mean(pf.xp(indpf,1),/double)*82.03125d0
print,ppf.POSITION_CENTER_OF_MASS_HALO(2,numhalof),mean(pf.xp(indpf,2),/double)*82.03125d0
print,ppf.VELOCITY_CENTER_OF_MASS_HALO(0,numhalof),mean(pf.vp(indpf,0),/double)*3.17682715523738E26/3.4995334104160499E17/100000.
print,ppf.VELOCITY_CENTER_OF_MASS_HALO(1,numhalof),mean(pf.vp(indpf,1),/double)*3.17682715523738E26/3.4995334104160499E17/100000.
print,ppf.VELOCITY_CENTER_OF_MASS_HALO(2,numhalof),mean(pf.vp(indpf,2),/double)*3.17682715523738E26/3.4995334104160499E17/100000.

xc=mean(pf.xp(indpf,0),/double)
yc=mean(pf.xp(indpf,1),/double)
zc=mean(pf.xp(indpf,2),/double)
print,ppf.RADIUS_MAXIMUM_HALO(numhalof)*ppf.NORMALISATION_RADIUS_PHYSICAL_HALO(numhalof),max(sqrt((pf.xp(indpf,0)-xc)^2+(pf.xp(indpf,1)-yc)^2+(pf.xp(indpf,2)-zc)^2))*3.17682715523738E26/3.08568E21*0.72
vxc=mean(pf.vp(indpf,0),/double)
vyc=mean(pf.vp(indpf,1),/double)
vzc=mean(pf.vp(indpf,2),/double)
print,ppf.VELOCITY_MAXIMUM_HALO(numhalof)*ppf.NORMALISATION_VELOCITY_HALO(numhalof),max(sqrt((pf.vp(indpf,0)-vxc+0.853085205526405*(pf.xp(indpf,0)-xc))^2  +(pf.vp(indpf,1)-vyc+0.853085205526405*(pf.xp(indpf,1)-yc))^2+(pf.vp(indpf,2)-vzc+0.853085205526405*(pf.xp(indpf,2)-zc))^2))*3.17682715523738d26/3.4995334104160499d17/100000.d0

print,ppf.DISPERSION_POSITION_HALO(numhalof)*ppf.NORMALISATION_RADIUS_PHYSICAL_HALO(numhalof),sqrt(mean((pf.xp(indpf,0)-xc)^2+(pf.xp(indpf,1)-yc)^2+(pf.xp(indpf,2)-zc)^2,/double))         *3.17682715523738E26/3.08568E21*0.72

print,ppf.DISPERSION_VELOCITY_HALO(numhalof)*ppf.NORMALISATION_VELOCITY_HALO(numhalof),sqrt(mean((pf.vp(indpf,0)-vxc+0.853085205526405*(pf.xp(indpf,0)-xc))^2  +(pf.vp(indpf,1)-vyc+0.853085205526405*(pf.xp(indpf,1)-yc))^2+(pf.vp(indpf,2)-vzc+0.853085205526405*(pf.xp(indpf,2)-zc))^2,/double))*3.17682715523738d26/3.4995334104160499d17/100000.d0

print,ppf.ENERGY_KINETIC_HALO(numhalof)*(0.5*ppf.NORMALISATION_VELOCITY_HALO(numhalof)^2*ppf.NORMALISATION_MASS_HALO(numhalof))/ppf.mass_halo(numhalof),0.5*mean((pf.vp(indpf,0)-vxc+0.853085205526405*(pf.xp(indpf,0)-xc))^2  +(pf.vp(indpf,1)-vyc+0.853085205526405*(pf.xp(indpf,1)-yc))^2+(pf.vp(indpf,2)-vzc+0.853085205526405*(pf.xp(indpf,2)-zc))^2,/double)*(3.17682715523738d26/3.4995334104160499d17/100000.d0)^2


epot=0.d0
eps=1.d0/128.d0/2.d0^6
epotmin=0.d0
imin=0L
epotitab=dblarr(nokf)
for i=0L,nokf-1 do begin
   d=sqrt((pf.xp(indpf(i),0)-pf.xp(indpf,0))^2+(pf.xp(indpf(i),1)-pf.xp(indpf,1))^2+(pf.xp(indpf(i),2)-pf.xp(indpf,2))^2+eps^2)
   unsurd=1.d0/d
   if i gt 0 then      epotinf=total(unsurd(0:i-1),/double) else epotinf=0.d0
   if i lt nokf-1 then epotsup=total(unsurd(i+1:nokf-1),/double) else epotsup=0.d0
   epoti=epotinf+epotsup
   epotitab(i)=epoti
   epot=epot+epoti
   if epoti ge epotmin then begin
      epotmin=epoti
      imin=i
   endif
   if i mod 1000 eq 0 then begin
      print,i,'/',nokf-1
      print,'epot=',epot
   endif
endfor

G=6.673099761655976d-8
mpart=1.d0/128.d0^3
unit_d = 3.39840378010395d-30
unit_l = 3.17682715523738d26
unit_t = 3.4995334104160499d17
epot=-0.5d0*epot*G*mpart*(unit_d*unit_l^3)/unit_l/nokf

print,ppf.ENERGY_SELF_BINDING_HALO(numhalof)*(0.5*ppf.NORMALISATION_VELOCITY_HALO(numhalof)^2*ppf.NORMALISATION_MASS_HALO(numhalof))/ppf.mass_halo(numhalof),epot*(1./100000.d0)^2
print,ppf.POSITION_MOST_BOUNDED_HALO(0,numhalof),pf.xp(indpf(imin),0)*82.03125d0
print,ppf.POSITION_MOST_BOUNDED_HALO(1,numhalof),pf.xp(indpf(imin),1)*82.03125d0
print,ppf.POSITION_MOST_BOUNDED_HALO(2,numhalof),pf.xp(indpf(imin),2)*82.03125d0

xx=pf.xp(indpf,0)-xc
yy=pf.xp(indpf,1)-yc
zz=pf.xp(indpf,2)-zc
vxx=pf.vp(indpf,0)-vxc
vyy=pf.vp(indpf,1)-vyc
vzz=pf.vp(indpf,2)-vzc

angx=yy*vzz-zz*vyy
angy=zz*vxx-xx*vzz
angz=xx*vyy-yy*vxx
print,ppf.ANGULAR_MOMENTUM_LAMBDA_PRIME_HALO(0,numhalof),total(angx,/double) /nokf /sqrt(2.d0) /(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.) /(sqrt(6.673099761655976E-8*(nokf/128.^3*3.39840378010395d-30*3.17682715523738d26^3)/(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.)/3.17682715523738E26))/(unit_t/unit_l)
print,ppf.ANGULAR_MOMENTUM_LAMBDA_PRIME_HALO(1,numhalof),total(angy,/double) /nokf /sqrt(2.d0) /(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.) /(sqrt(6.673099761655976E-8*(nokf/128.^3*3.39840378010395d-30*3.17682715523738d26^3)/(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.)/3.17682715523738E26))/(unit_t/unit_l)
print,ppf.ANGULAR_MOMENTUM_LAMBDA_PRIME_HALO(2,numhalof),total(angz,/double) /nokf /sqrt(2.d0) /(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.) /(sqrt(6.673099761655976E-8*(nokf/128.^3*3.39840378010395d-30*3.17682715523738d26^3)/(nokf/128.^3/(4./3.*!pi*1.*200.))^(1./3.)/3.17682715523738E26))/(unit_t/unit_l)

print,ppf.ANGULAR_MOMENTUM_LAMBDA_HALO(0,numhalof)/ppf.ANGULAR_MOMENTUM_LAMBDA_PRIME_HALO(0,numhalof),sqrt(abs(ppf.ENERGY_SELF_BINDING_HALO(numhalof)+ppf.ENERGY_KINETIC_HALO(numhalof)))

print,sqrt(total(ppf.INERTIA_EIGEN_VALUES_HALO(*,numhalof)^2,/double)),ppf.DISPERSION_POSITION_HALO(numhalof) ;just check consistency

;rem do not check eigenvector nor cos angle-> too tired. hopefuly
;ok. I tested long time ago I think...


;CHECK PROF
xcpot=pf.xp(indpf(imin),0)
ycpot=pf.xp(indpf(imin),1)
zcpot=pf.xp(indpf(imin),2)

delta_lr = 0.1d0
rho=dblarr(pff.nprof)
rho_integrated=dblarr(pff.nprof)
minfr=dblarr(pff.nprof)
r=double(pff.PROFILE_RADIAL_BINS_HALO/0.72)
rmin=exp(alog(pff.PROFILE_RADIAL_BINS_HALO/0.72)-delta_lr/2.)
rmax=exp(alog(pff.PROFILE_RADIAL_BINS_HALO/0.72)+delta_lr/2.)
r=r/(82.03125d0*1000./0.72*0.903640895892252)
rmin=rmin/(82.03125d0*1000./0.72*0.903640895892252)
rmax=rmax/(82.03125d0*1000./0.72*0.903640895892252)

d=sqrt((pf.xp(indpf,0)-xcpot)^2+(pf.xp(indpf,1)-ycpot)^2+(pf.xp(indpf,2)-zcpot)^2)

ind=where(d lt rmin(0),nok)
nok0=nok
for i=0L,pff.nprof -1 do begin
   ind=where(d ge rmin(i) and d lt rmax(i),nok)
   rho(i)=nok
   if i gt 0 then  minfr(i)=minfr(i-1)+nok else minfr(i)=nok0+nok
   rho_integrated(i)= minfr(i)*mpart/(4./3.*!pi*rmax(i)^3)
endfor
print,total(rho,/double)+nok0,minfr(pff.nprof -1),nokf
mpart=1.d0/128.d0^3
rho=rho*mpart/(4.d0/3.d0*!dpi*(rmax^3-rmin^3))
tek_color
plot,pff.PROFILE_RADIAL_BINS_HALO/0.72,pff.PROFILE_DENSITY_HALO(*,numhalof),/xlog,/ylog,yr=[1d-2,1d6],psym=-4
oplot,pff.PROFILE_RADIAL_BINS_HALO/0.72,rho,psym=-4, color=3,lines=2

stop

plot,pff.PROFILE_RADIAL_BINS_HALO/0.72,pff.PROFILE_DENSITY_INTEGRATED_HALO(*,numhalof),/xlog,/ylog,yr=[1d-2,1d7],psym=-4
oplot,pff.PROFILE_RADIAL_BINS_HALO/0.72,rho_integrated,psym=-4, color=3,lines=2

stop
plot,pff.PROFILE_RADIAL_BINS_HALO/0.72,pff.PROFILE_CIRCULAR_VELOCITY_HALO(*,numhalof),/xlog,/ylog,psym=-4
vcirc=(sqrt(6.673099761655976E-8*(minfr/128.^3*3.39840378010395d-30*3.17682715523738d26^3)/rmax/3.17682715523738E26))/100000.
oplot,pff.PROFILE_RADIAL_BINS_HALO/0.72,vcirc,psym=-4, color=3,lines=2
ind=where(vcirc eq max(vcirc))
print,minfr(ind(0))*mpart*3.39840378010395d-30*3.17682715523738d26^3/1.98892E33*0.72,rmax(ind(0))*(82.03125d0*1000./0.72*0.903640895892252)*0.72,vcirc(ind(0)),rho_integrated(ind(0))
print,ppf.CIRCULAR_VALUES_HALO(*,numhalof)

ind=where(abs(rho_integrated-200.) eq min(abs(rho_integrated-200.)))
print,ppf.delta_values_halo(0,numhalof),ppf.delta_values_halo(1,numhalof)/0.72
;print,rho_integrated(ind(0))*4.d0/3.d0*!dpi*rmax(ind(0))^3*unit_d*unit_l^3*0.72d0/1.98892E33, pff.PROFILE_RADIAL_BINS_HALO(ind(0))/0.72*exp(delta_lr/2.)
rm=exp(interpol(alog(rmax),alog(rho_integrated),alog(200.)))
print,200.*4.d0/3.d0*!dpi*rm^3*unit_d*unit_l^3*0.72d0/1.98892E33,rm*(82.03125d0*1000./0.72*0.903640895892252)
stop


delta_lr = 0.2d0
rho=dblarr(pff.nprof178)
rho_integrated=dblarr(pff.nprof178)
minfr=dblarr(pff.nprof178)
r=double(pff.PROFILE_RADIAL_BINS_RDELTA_HALO)*ppf.delta_values_halo(1,numhalof)/0.72d0
rmin=exp(alog(double(pff.PROFILE_RADIAL_BINS_RDELTA_HALO)*double(ppf.delta_values_halo(1,numhalof))/0.72d0)-delta_lr/2.d0)
rmax=exp(alog(double(pff.PROFILE_RADIAL_BINS_RDELTA_HALO)*double(ppf.delta_values_halo(1,numhalof))/0.72d0)+delta_lr/2.d0)
r=r/(82.03125d0*1000.d0/0.72d0*0.903640895892252d0)
rmin=rmin/(82.03125d0*1000.d0/0.72d0*0.903640895892252d0)
rmax=rmax/(82.03125d0*1000.d0/0.72d0*0.903640895892252d0)

;d=sqrt((pf.xp(indpf,0)-xcpot)^2+(pf.xp(indpf,1)-ycpot)^2+(pf.xp(indpf,2)-zcpot)^2)
d=sqrt((pf.xp(indpf,0)-ppf.POSITION_MOST_BOUNDED_HALO(0,numhalof)/82.03125d0)^2+(pf.xp(indpf,1)-ppf.POSITION_MOST_BOUNDED_HALO(1,numhalof)/82.03125d0)^2+(pf.xp(indpf,2)-ppf.POSITION_MOST_BOUNDED_HALO(2,numhalof)/82.03125d0)^2)

ind=where(d lt rmin(0),nok)
nok0=nok
for i=0L,pff.nprof178 -1 do begin
   ind=where(d ge rmin(i) and d lt rmax(i),nok)
   rho(i)=nok
   if i gt 0 then  minfr(i)=minfr(i-1)+nok else minfr(i)=nok0+nok
   rho_integrated(i)= minfr(i)*mpart/(4.d0/3.d0*!dpi*rmax(i)^3)
endfor
print,total(rho,/double)+nok0,minfr(pff.nprof178 -1),nokf
mpart=1.d0/128.d0^3
rho=rho*mpart/(4.d0/3.d0*!dpi*(rmax^3-rmin^3))
tek_color
plot,pff.PROFILE_RADIAL_BINS_RDELTA_HALO,pff.PROFILE_DENSITY_BINS_RDELTA_HALO(*,numhalof),/xlog,/ylog,yr=[1d-2,1d6],psym=-4
oplot,pff.PROFILE_RADIAL_BINS_RDELTA_HALO,rho,psym=-4, color=3,lines=2
stop
plot,pff.PROFILE_RADIAL_BINS_RDELTA_HALO,pff.PROFILE_DENSITY_INTEGRATED_BINS_RDELTA_HALO(*,numhalof),/xlog,/ylog,yr=[1d-2,1d7],psym=-4
oplot,pff.PROFILE_RADIAL_BINS_RDELTA_HALO,rho_integrated,psym=-4, color=3,lines=2



print,'END'
stop
end

pro test_sod_vs_haloanalyzer
rd_pfof_halo_snaporcone_part_prop_hdf5,pps,pfs,file='/obs/yrasera/proj/fulluniverse/code/haloanalyzer_new/psod_halo_snap_part_prop_boxlen82.03125_n128_lcdmw7.h5'

for i=0L,pps.nhalo-1L do begin
   plot,pfs.PROFILE_RADIAL_BINS_RDELTA_HALO*exp(0.1),pfs.PROFILE_DENSITY_INTEGRATED_BINS_RDELTA_HALO(*,i),/xlog,/ylog,psym=-4,xr=[0.01,10],yr=[1,1e6]
   oplot,[1,1],[1d-10,1d10]
   oplot,[1d-10,1d10],[200,200]
   wait,0.1
   if i mod 100 eq 0 then print,'i=',i,'/',pps.nhalo-1L
end

;for i=0L,pps.nhalo-1L do begin
;   plot,pfs.PROFILE_RADIAL_BINS_RDELTA_HALO,pfs.PROFILE_DENSITY_BINS_RDELTA_HALO(*,i),/xlog,/ylog,psym=-4,xr=[0.01,10],yr=[1,1e6]
;   oplot,[1,1],[1d-10,1d10]
;   oplot,[1d-10,1d10],[200,200]
;   stop
;end

end



pro test_fof_vs_haloanalyzer
rd_pfof_halo_snaporcone_part_prop_hdf5,ppf,pff,file='/obs/yrasera/proj/fulluniverse/code/haloanalyzer_new/pfof_halo_snap_part_prop_boxlen82.03125_n128_lcdmw7.h5'

for i=0L,ppf.nhalo-1L do begin
   plot,pff.PROFILE_RADIAL_BINS_RDELTA_HALO*exp(0.1),pff.PROFILE_DENSITY_INTEGRATED_BINS_RDELTA_HALO(*,i),/xlog,/ylog,psym=-4,xr=[0.01,10],yr=[1,1e6]
   oplot,[1,1],[1d-10,1d10]
   oplot,[1d-10,1d10],[200,200]
   wait,0.1
   if i mod 100 eq 0 then print,'i=',i,'/',ppf.nhalo-1L
end

;for i=0L,ppf.nhalo-1L do begin
;   plot,pff.PROFILE_RADIAL_BINS_RDELTA_HALO,pff.PROFILE_DENSITY_BINS_RDELTA_HALO(*,i),/xlog,/ylog,psym=-4,xr=[0.01,10],yr=[1,1e6]
;   oplot,[1,1],[1d-10,1d10]
;   oplot,[1d-10,1d10],[200,200]
;   stop
;end

end
