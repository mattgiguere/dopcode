pro temp_corl, acena=acena, acenb=acenb, track=track,guide=guide

; beginning in Mar 2010, air temp was monitored in the Coude room at
; CTIO and recorded in the FITS headers. 

; this program will check for correlations between detrended (linear
; trend removed) velocities in alpha Cen A and B

; RESTORE VST for aCEN A and remove linear trend
restore,'/home/fischer/dop2/vank/vstbank/vst128620.dat'  & cfa=cf3
x=where(cfa.jd gt 15200.)   ; air temp not in headers before this JD
cfa=cfa[x]  ;post temp monitor installation
coef1=polyfit(cfa.jd,cfa.mnvel,2)
plot,cfa.jd,cfa.mnvel,ps=8
oplot,cfa.jd,poly(cfa.jd,coef1)
cfa.mnvel=cfa.mnvel-poly(cfa.jd,coef1)
cfa.mnvel=cfa.mnvel-median(cfa.mnvel)   ; recenter velocities about zero

; RESTORE VST for aCEN B and remove linear trend
restore,'/home/fischer/dop2/vank/vstbank/vst128621.dat'  & cfb=cf3
x=where(cfb.jd gt 15200.)   ; air temp not in headers before this JD
cfb=cfb[x]
coef2=polyfit(cfb.jd,cfb.mnvel,2)
plot,cfb.jd,cfb.mnvel,ps=8
oplot,cfb.jd,poly(cfb.jd,coef2)
cfb.mnvel=cfb.mnvel-poly(cfb.jd,coef2)
cfb.mnvel=cfb.mnvel-median(cfb.mnvel)   ; recenter velocities about zero

!p.charsize=1.8

; GUIDING TEST: images below were guided center, left and right on May 8, 2010
if keyword_set(guide) then begin
   restore,'vstbank/vst128620.dat'
   cobnm1='rqa15.'+strcompress(string(indgen(30)+5448),/rem)
   cobnm2='rqa15.'+strcompress(string(indgen(30)+5553),/rem)
   lobnm1='rqa15.'+strcompress(string(indgen(30)+5483),/rem)
   lobnm2='rqa15.'+strcompress(string(indgen(30)+5588),/rem)
   robnm1='rqa15.'+strcompress(string(indgen(30)+5518),/rem)
   robnm2='rqa15.'+strcompress(string(indgen(30)+5623),/rem)

   cobnm=[cobnm1,cobnm2]  & cnum=n_elements(cobnm)
   lobnm=[lobnm1,lobnm2]  & lnum=n_elements(lobnm)
   robnm=[robnm1,robnm2]  & rnum=n_elements(robnm)

   for i=0,cnum-1 do begin
      x=where(cobnm[i] eq cf3.obnm,nx)
      if i eq 0 then ccf3=cf3[x]
      if i gt 0 and nx gt 0 then ccf3=[ccf3,cf3[x]]
   endfor

   for i=0,lnum-1 do begin
      x=where(lobnm[i] eq cf3.obnm,nx)
      if i eq 0 then lcf3=cf3[x]
      if i gt 0 and nx gt 0 then lcf3=[lcf3,cf3[x]]
   endfor

   for i=0,rnum-1 do begin
      x=where(robnm[i] eq cf3.obnm,nx)
      if i eq 0 then rcf3=cf3[x]
      if i gt 0 and nx gt 0 then rcf3=[rcf3,cf3[x]]
   endfor

   plot,ccf3.jd-15325.54,ccf3.mnvel,ps=8,col=1, xra=[0.18,0.35],/xsty,yra=[70,140],$
        /ysty,titl='Guiding test, aCen A', ytit='RV [m s!u -1 !n]',xtit='!6 Time [days]
   oplot,lcf3.jd-15325.54,lcf3.mnvel,ps=8,col=90
   oplot,rcf3.jd-15325.54,rcf3.mnvel,ps=8,col=222
   xyouts,0.19,130,/data,'!6 RMS (center): '+strmid(strcompress(string(stddev(ccf3.mnvel)),/rem),0,4)
   !p.color=90
   xyouts,0.19,125,/data,'!6 RMS (left): '+strmid(strcompress(string(stddev(lcf3.mnvel)),/rem),0,4)
   !p.color=222
   xyouts,0.19,120,/data,'!6 RMS (right): '+strmid(strcompress(string(stddev(rcf3.mnvel)),/rem),0,4)
   print,'Mean chi (center): ',mean(ccf3.mdchi),' Mean err (center): ',mean(ccf3.errvel)
   print,'Vel RMS (center): ',stddev(ccf3.mnvel)
   print,'Mean chi (left): ',mean(lcf3.mdchi),' Mean err (left): ',mean(lcf3.errvel)
   print,'Vel RMS (left): ',stddev(lcf3.mnvel)
   print,'Mean chi (right): ',mean(rcf3.mdchi),' Mean err (right): ',mean(rcf3.errvel)
   print,'Vel RMS (right): ',stddev(rcf3.mnvel)
 
   stop
endif   ;GUIDING test


; Consider aCenA and B seperately or (if acena or acenb not set) as merged sets
if keyword_set(acena) then cf=cfa
if keyword_set(acenb) then cf=cfb
if ~keyword_set(acena) and ~keyword_set(acenb) then cf=[cfa,cfb]
i=sort(cf.jd)
cf=cf[i]

; Keyword TRACK
; DO mean RV's for aCENA and aCENB track each other from night to night?
if keyword_set(track) then begin
   jd0a=fix(cfa[0].jd)
   velplot,cfa,9./24.,d,v,e,n,c,duma  &  cfa=duma
   velplot,cfb,9./24.,d,v,e,n,c,dumb  &  cfb=dumb
!p.charsize=2.4
!p.font=1
   plot,duma.jd-jd0a,duma.mnvel,ps=4,col=1,yra=[-20,20],xtit='!6 JD - 2455200.',ytit='!6 RV m s!u-1!n'
   oplot,dumb.jd-jd0a,dumb.mnvel,ps=8,col=222
stop
endif

;stop
num=n_elements(cf)  ;number of observations
tair=fltarr(num)
ctemp=fltarr(num)
ha_dec=fltarr(num)
vel=fltarr(num)
err=fltarr(num)
chi=fltarr(num)
hour=intarr(num) 
cts=lonarr(num)
par0=fltarr(num)  
par1=fltarr(num)  &  par2=fltarr(num)   & par3=fltarr(num)  & par4=fltarr(num) 
par5=fltarr(num)  &  par6=fltarr(num)   & par7=fltarr(num)  & par8=fltarr(num) 


if ~keyword_set(acena) and ~keyword_set(acenb) then begin
   print,'if looking for correlations with TAIR, specifiy acena or acenb'
   stop
endif

for i=0,num-1 do begin
   filenm=cf[i].obnm+'.fits'
   hd=headfits('/mir7/fitspec/'+filenm)
      tair[i]=sxpar(hd,'TEMPAIR1')
      ctemp[i]=sxpar(hd,'CCDTEMP')
      ha=sxpar(hd,'HA')
      if strmid(ha,0,1) eq '-' then hour[i]=strmid(ha,0,3)
      if strmid(ha,0,1) ne '-' then hour[i]=strmid(ha,0,2)
      hour[i]=fix(hour[i])
      min=strmid(ha,4,2)  &  min=fix(min)
      sec=strmid(ha,7,2)  &  sec=float(sec)
      ha_dec[i]=ten(hour[i],min,sec)
      vel[i]=cf[i].mnvel
      err[i]=cf[i].errvel
      chi[i]=cf[i].mdchi
      cts[i]=cf[i].cts
   restore,'/mir7/files/vdb*'+cf[i].obnm
   par0[i]=mean(vd[180:182].psf[60])
   par1[i]=mean(vd[180:182].psf[45])
   par2[i]=mean(vd[180:182].psf[50])
   par3[i]=mean(vd[180:182].psf[55])
   par4[i]=mean(vd[180:182].psf[58])
   par5[i]=mean(vd[180:182].psf[62])
   par6[i]=mean(vd[180:182].psf[65])
   par7[i]=mean(vd[180:182].psf[70])
   par8[i]=mean(vd[180:182].psf[75])
endfor

; vel as a function of HA
x=where(ha_dec lt 0.)
plot,ha_dec,vel,ps=8,xtit='HA',ytit='Vel [m/s]',/xsty,/ysty  
coef=poly_fit(ha_dec,vel,1)  &  oplot,ha_dec,poly(ha_dec,coef)
stop

plot,ha_dec[x],vel[x],ps=8,xtit='HA',ytit='Vel [m/s]',/xsty,/ysty  
coef1=poly_fit(ha_dec[x],vel[x],1)  &  oplot,ha_dec[x],poly(ha_dec[x],coef1)
stop

plot,err,vel,ps=8,xtit='Errvel',ytit='Vel [m/s]',/xsty,/ysty  
coef2=poly_fit(err,vel,1)  &  oplot,err,poly(err,coef2)
stop

plot,chi,vel,ps=8,xtit='Chi',ytit='Vel [m/s]',/xsty,/ysty  
coef3=poly_fit(chi,vel,1)  &  oplot,chi,poly(chi,coef3)
stop

plot,tair,vel,ps=8,xtit='TAIR [C]',ytit='Vel [m/s]',xra=[10,20],/xsty,/ysty 
coef4=poly_fit(tair,vel,1)  &  oplot,tair,poly(tair,coef4)
stop

plot,ctemp,vel,ps=8,xtit='CCD TEMP [C]',ytit='Vel [m/s]',/xsty,/ysty 
coef5=poly_fit(ctemp,vel,1)  &  oplot,ctemp,poly(ctemp,coef5)
stop

plot,sqrt(cts),vel,ps=8,xtit='SNR',ytit='Vel [m/s]',/xsty,/ysty 
coef6=poly_fit(sqrt(cts),vel,1)  &  oplot,sqrt(cts),poly(sqrt(cts),coef6)
stop

plot,ha_dec[x],sqrt(cts),ps=8,xtit='HA',ytit='SNR',/xsty,/ysty
stop

; err
plot,sqrt(cts),err,ps=8,xtit='SNR',ytit='Err [m/s]',/xsty,/ysty 
coef7=poly_fit(sqrt(cts),err,1)  &  oplot,sqrt(cts),poly(sqrt(cts),coef7)
stop

!p.multi=[0,3,3]
plot,par1,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 1'
   i1=sort(par1)  &  vel1=vel[i1]   &   par1=par0[i1]
   coef1=poly_fit(par1,vel1,1)
   xarr1=findgen(10)/(max(par1)-min(par1)) + min(par1)
   oplot,xarr1,poly(xarr1,coef1)
plot,par2,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 2'
   i2=sort(par2)  &  vel2=vel[i2]   &   par2=par2[i2]
   coef2=poly_fit(par2,vel2,1)
   xarr2=findgen(10)/(max(par2)-min(par2)) + min(par2)
   oplot,xarr2,poly(xarr2,coef2)
plot,par3,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 3'
   i3=sort(par3)  &  vel3=vel[i3]   &   par3=par3[i3]
   coef3=poly_fit(par3,vel3,1)
   xarr3=findgen(10)/(max(par3)-min(par3)) + min(par3)
   oplot,xarr3,poly(xarr3,coef3)
plot,par4,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 4'
   i4=sort(par4)  &  vel4=vel[i4]   &   par4=par4[i4]
   coef4=poly_fit(par4,vel4,1)
   xarr4=findgen(10)/(max(par4)-min(par4)) + min(par4)
   oplot,xarr4,poly(xarr4,coef4)
plot,par0,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 0'
   i0=sort(par0)  &  vel0=vel[i0]   &   par0=par0[i0]
   coef0=poly_fit(par0,vel0,1)
   xarr0=findgen(10)/(max(par0)-min(par0)) + min(par0)
   oplot,xarr0,poly(xarr0,coef0)
plot,par5,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 5'
   i5=sort(par5)  &  vel5=vel[i5]   &   par5=par5[i5]
   coef5=poly_fit(par5,vel5,1)
   xarr5=findgen(10)/(max(par5)-min(par5)) + min(par5)
   oplot,xarr5,poly(xarr5,coef5)
plot,par6,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 6'
   i6=sort(par6)  &  vel6=vel[i6]   &   par6=par6[i6]
   coef6=poly_fit(par6,vel6,1)
   xarr6=findgen(10)/(max(par6)-min(par6)) + min(par6)
   oplot,xarr6,poly(xarr6,coef6)
plot,par7,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 7'
   i7=sort(par7)  &  vel7=vel[i7]   &   par7=par7[i7]
   coef7=poly_fit(par7,vel7,1)
   xarr7=findgen(10)/(max(par7)-min(par7)) + min(par7)
   oplot,xarr7,poly(xarr7,coef7)
plot,par8,vel,ps=8,yra=[-40,40],/ysty,/xsty,symsize=0.4,titl='!6 PAR 8'
   i8=sort(par8)  &  vel8=vel[i8]   &   par8=par8[i8]
   coef8=poly_fit(par8,vel8,1)
   xarr8=findgen(10)/(max(par8)-min(par8)) + min(par8)
   oplot,xarr8,poly(xarr8,coef8)

stop

end
