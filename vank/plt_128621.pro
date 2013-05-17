pro plt_128621

restore,'/home/fischer/dop2/vank/vst128621_gd.dat'
cf3=cf3_gd

velplot,cf3,15./24.,d,v,e,n,c,dum
coef=poly_fit(dum.jd,dum.mnvel,1)
cf3.mnvel=cf3.mnvel-poly(cf3.jd,coef)

;plothist,cf3.mnvel
;x=where(abs(cf3.mnvel) le 18)
;cf3=cf3[x]
xx=where(cf3.errvel lt 1.5*median(cf3.errvel))
cf3=cf3[xx] 
velplot,cf3,titl='!6 Alpha Cen B'


openw,1,'acenb.txt'
for i=0,n_elements(cf3.mnvel)-1 do printf,1,cf3[i].jd,cf3[i].mnvel,cf3[i].errvel,f='(f16.5,f9.1,f9.1)'
close,1

stop
velplot,cf3,15./24.,d,v,e,c,n,dum
velplot,dum 
stop

end
