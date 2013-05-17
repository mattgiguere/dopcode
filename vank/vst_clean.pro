pro vst_clean, star=star, alerr=alerr

if ~keyword_set(alerr) then alerr=2.  ; allowed error range

if star eq '128621' then vstnm='vst128621.dat'
if star eq '128620' then vstnm='vst128620.dat'

restore,'vstbank/'+vstnm
flag=[1]

minjd=min(cf3.jd) & minjd=fix(minjd)   &  minjd=minjd*1.0
maxjd=max(cf3.jd) & maxjd=fix(maxjd)   &  maxjd=maxjd*1.0
range=fix(maxjd-minjd) 

for i=0,range-1 do begin
   x=where(fix(cf3.jd) eq minjd+i,nx)
   if nx gt 0 then begin
      tmpflag=intarr(nx)
      mderr=median(cf3[x].errvel)
      mdvel=median(cf3[x].mnvel)
      xx=where(cf3[x].mnvel gt mdvel-(alerr*mderr) and cf3[x].mnvel lt mdvel+(alerr*mderr))
      if i eq 0 then plothist,cf3[x[xx]].mnvel-mdvel,xra=[-80,80],yra=[0,50],bin=3
      if i gt 0 then plothist,cf3[x[xx]].mnvel-mdvel,/overplot, col=i*7,bin=3
      tmpflag[xx]=1
      flag=[flag,tmpflag]
   endif
endfor
;stop
num=n_elements(flag)
flag=flag[1:num-1] ; remove the first dummy index 
cf3_gd=cf3[where(flag eq 1)]
velplot,cf3_gd 
save,cf3_gd,f='vst'+star+'_gd.dat'

stop
end
