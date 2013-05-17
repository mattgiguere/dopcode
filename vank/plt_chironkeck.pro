pro plt_chironkeck, hdcopy=hdcopy

restore,'/tous/mir3/vel/vst10700.dat'  ;cf3 for keck
x=where(cf3.mnvel lt 10 and cf3.mnvel gt -10)  & cf3=cf3[x]
cfk=cf3[x]
y=where(cfk.dewar eq 102)
z=where(cfk.dewar eq 103)
print,'stdev new keck: ',stddev(cfk[z].mnvel)

restore,'vst10700_nov2011.dat'
velplot,cf3,4./24.,d,v,e,n,c,cfo

restore,'vstbank/vst10700.dat'
velplot,cf3,0.0001,d,v,e,n,c,cfn  
velplot,cfn,4./24.,d,v,e,n,c,cfc  
print,'stddev CHIRON: ',stddev(cfc.mnvel)

!p.font=1
!x.charsize=1.9
!y.charsize=1.7

if keyword_set(hdcopy) then ps_open,'chiron_vs_keck',/color
plot,sqrt(cfc.cts),cfc.errvel,ps=8,yra=[0,3],xra=[50,550],/xsty,$
     xtitl='!6 SNR',ytitl='!6 Errvel [m s!u-1!n]',/nodata
oplot,sqrt(cfc.cts),cfc.errvel,ps=6,thick=4  
oplot,sqrt(cfk[y].cts),cfk[y].errvel,ps=5,col=160 
oplot,sqrt(cfk[z].cts),cfk[z].errvel,ps=5,col=90  
oplot,sqrt(cfo.cts), cfo.errvel,ps=7,col=210
plots,[150,450],[1.0,1.0],linesty=2  
;plots,[150,450],[0.5,0.5],linesty=2,col=222,thick=2

plots,[350,350],[2.65,2.65],ps=5,col=160,thick=6,symsize=2  
xyouts,360,2.6,'!6 Keck pre-upgrade',/data,size=1.9

plots,[350,350],[2.45,2.45],ps=5,col=90,thick=6,symsize=2  
xyouts,360,2.4,'!6 Keck post-upgrade',/data,size=1.9

plots,[350,350],[2.25,2.25],ps=6,col=0,thick=4,symsize=2  
xyouts,360,2.2,'!6 CHIRON',/data,size=1.9

cts=cfc.cts
snr=sqrt(cts)
ii=sort(snr)
y=cfc.errvel
snr=snr[ii]
yrr=y[ii]

; snr=[130.046, 176.156, 212.770, 225.424, 262.025, 270.104, 294.297, 374.015, 390.377, 414.402,$
;      422.882, 432.744, 439.030, 444.672, 520.773]
; yrr=[0.977008, 0.880636, 0.749620, 0.583097, 0.450807, 0.453540, 0.775650, 0.378712, 0.355396, 0.330641,$
;     0.368585, 0.308208, 0.384679, 0.342034, 0.270267]

 err=fltarr(n_elements(snr))+0.05
 err[3]=4.
 err[8]=3.5

pinit=[0.5,-0.006]
lam=mpfitfun('expdecay',snr,yrr,err,pinit)
xsnr=findgen(50)*10+100.
ny=expdecay(xsnr,lam)
oplot,xsnr,ny,col=222,thick=3,linesty=2


if keyword_set(hdcopy) then ps_close


g1=where(cfk.jd gt 16110)
g2=where(cfc.jd gt 16110)
velplot,cfk[g1],4./24.,dk,vk,ek,nk,cj,cfkc
velplot,cfc[g2],4./24.,dc,vc,ec,nc,cc,cfcc

plot,cfkc.jd,cfkc.mnvel+4.,yra=[-20,20],ps=5,/nodata
oplot,cfkc.jd,cfkc.mnvel+4.,ps=5,col=90
oplot,cfcc.jd,cfcc.mnvel,ps=6

stop
end
