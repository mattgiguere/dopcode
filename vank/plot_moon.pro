pro plot_moon

restore,'vstbank/vstmoon.dat'
num=n_elements(cf3)

diff=fltarr(num)
for i=0,num-2 do diff[i] = cf3[i+1].jd - cf3[i].jd
break = where(diff gt 0.1,nx) 
;stop

for i=0,nx-1 do begin
   if i eq 0 then begin
      cf=cf3[0:break[0]]
      jd0=14929.-cf[0].jd
      cf.jd=cf.jd+jd0
   endif
   if i gt 0 then begin
      cf = cf3[break[i-1]+1:break[i]] 
;      jd0=cf[0].jd
   endif
;   cf.jd = cf.jd-jd0
   coef=poly_fit(cf.jd,cf.mnvel,1)
   cf.mnvel=cf.mnvel-poly(cf.jd,coef)
   velplot,cf, 0.167/24.  ; 10 minute bins
   if i eq 0 then cf_all = cf else cf_all = [cf, cf_all]
;   stop
endfor

i=sort(cf_all.jd)
cf_all=cf_all[i]
velplot,cf_all, 0.167/24. 

stop
end
