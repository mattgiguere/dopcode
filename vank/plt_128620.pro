pro plt_128620

restore,'/home/fischer/dop2/vank/vst128620_gd.dat'
cf3=cf3_gd

velplot,cf3,15./24.,d,v,e,n,c,dum
coef=poly_fit(dum.jd,dum.mnvel,2)
cf3.mnvel=cf3.mnvel-poly(cf3.jd,coef)

;x=where(abs(cf3.mnvel) lt 15)
;cf3=cf3[x]
xx=where(cf3.errvel lt 1.5*median(cf3.errvel))
cf3=cf3[xx] 
velplot,cf3


openw,1,'acena.txt'
for i=0,n_elements(cf3.mnvel)-1 do printf,1,cf3[i].jd,cf3[i].mnvel,cf3[i].errvel,f='(f16.5,f9.1,f9.1)'
close,1

stop
velplot,cf3,(15./60)/24.,d,v,e,n,c,dum
velplot,dum
stop

end
