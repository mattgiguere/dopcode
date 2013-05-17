pro psf_rv, star=star

; look for correlations between psf paramters and velocities
; dfischer jul 2009

if keyword_set(star) then restore,'vstbank/vst'+star+'.dat' else $
   stop,'enter the starname, (e.g., star=128621)'

nobs=n_elements(cf1)
clight=299792458.

   md_psfpar=fltarr(nobs,11)
   vel=fltarr(nobs)
   obnm=strarr(nobs) 
   nchunks=n_elements(cf1.obnm)
   namp=17
   

for i=0,nobs-1 do begin
print,i
   cfile='cdamb'+star+'_'+cf1[i].obnm
   obnm[i]=cf1[i].obnm
   restore,'/tous/mir7/files/'+cfile
   vel[i]=cf1[i].mnvel
   amp=fltarr(nchunks,namp) 
   for j=0,nchunks-1 do begin
      for k=0,namp-1 do begin
         amp[j,k]=(*chunk_arr[j].free_par)[1].amp[k] 
      endfor
   endfor
   for l=0,10 do md_psfpar[i,l]=median( amp[*,l] )
endfor

save,obnm,md_psfpar,vel,f='psf_vel_'+star+'.dat'

!p.font=1
!p.charsize=1.5
!x.charsize=1.1
!y.charsize=1.1
!p.multi=[0,2,1]
ps_open,'psf-'+star+'_1-10'
plot,md_psfpar[*,1],vel,ps=8,xtit='!6PSF par[1]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
plot,md_psfpar[*,10],vel,ps=8,xtit='!6PSF par[10]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
ps_close

ps_open,'psf-'+star+'_2-9'
plot,md_psfpar[*,2],vel,ps=8,xtit='!6PSF par[2]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
plot,md_psfpar[*,9],vel,ps=8,xtit='!6PSF par[9]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
ps_close

ps_open,'psf-'+star+'_3-8'
plot,md_psfpar[*,3],vel,ps=8,xtit='!6PSF par[3]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
plot,md_psfpar[*,8],vel,ps=8,xtit='!6PSF par[8]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
ps_close

ps_open,'psf-'+star+'_4-7'
plot,md_psfpar[*,4],vel,ps=8,xtit='!6PSF par[4]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
plot,md_psfpar[*,7],vel,ps=8,xtit='!6PSF par[7]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
ps_close

ps_open,'psf-'+star+'_5-6'
plot,md_psfpar[*,5],vel,ps=8,xtit='!6PSF par[5]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
plot,md_psfpar[*,6],vel,ps=8,xtit='!6PSF par[6]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
ps_close

ps_open,'psf-'+star+'_0'
plot,md_psfpar[*,0],vel,ps=8,xtit='!6PSF par[0]',ytit='!6Mean Vel [m s!u-1!n]',symsize=0.5,yra=[-40,40]
ps_close

!p.multi=[0,1,1]

stop

end

