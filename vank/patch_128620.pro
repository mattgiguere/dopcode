pro patch_128620

restore,'psf_vel_128620.dat'  ;obnm, md_psfpar, vel 
par0=md_psfpar[*,0]

i=sort(par0)
par0=par0[i]
vel=vel[i]
obnm=obnm[i]
num=n_elements(par0)

coef=poly_fit(par0,vel,2)
!p.multi=[0,2,1]
plot,par0,vel,ps=8
oplot,par0, poly(par0,coef),col=222,thick=5

; dvel = coef[0] + coef[1]*par0
; cf3.mnvel = cf3.mnvel - dvel 

dvel = (coef[0] + coef[1]*par0 + coef[2]*par0^2)
nvel = vel - dvel 
plot,par0, nvel,ps=8

!p.multi=[0,1,1]

restore,'vstbank/vst128620.dat'
for i=0, num-1 do begin
   x=where(obnm[i] eq cf1.obnm,nx)
   if nx eq 0 then stop
   cf3[x].mnvel = cf3[x].mnvel - dvel[i]
endfor

velplot,cf3
save,cf1,cf3,f='vstbank/vst128620.dat


stop
end
