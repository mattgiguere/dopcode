pro clean_acen, tag=tag, star=star

vdtag='vd'+tag
vank,star,tag,vdtag, cf1,cf3,/ctio, mct=50000.
restore,'vstbank/vst'+star+'.dat'
velplot,cf3

flag=intarr(n_elements(cf3))+1
bd=['achi120716.1125','achi120716.1126','achi120716.1128','achi120716.1129']
bd=[bd,'achi120808.0000','achi120808.0001','achi120808.0002','achi120808.0003','achi120808.0004']
for i=0,n_elements(cf3)-1 do begin
	x=where(cf3[i].obnm eq bd,nx)
	if nx eq 1 then flag[i]=0
endfor
g1=where(flag eq 1) 
cf3=cf3[g1]

thresh_RMS = 3.
xdate=fix(cf3.jd) 
g=uniq(xdate)
xdate=xdate[g]
if tag eq 'am' then begin
	xk=where(xdate ne 16131)  
	xdate=xdate[xk]
	endif
if tag eq 'c' then begin 
	xk=where(xdate ne 16126)  
	xdate=xdate[xk]
endif
	
num=n_elements(xdate) 

good=[0]
sdev=fltarr(num)+999.

for i=0,num-1 do begin
   x=where(cf3.jd gt xdate[i] and cf3.jd lt xdate[i]+1.,nx) 
   if nx gt 2 then begin
      if stddev(cf3[x].mnvel) lt thresh_RMS then begin
         good=[good,x]
         sdev[i]=stddev(cf3[x].mnvel)
         print,sdev[i]
      endif
   endif
endfor 
good=good[1:n_elements(good)-1]
xs=where(sdev ne 999.)  &  sdev=sdev[xs]

plothist,sdev,bin=0.2

;print,'Plotting pergram of data with slope included' 
cf=cf3[good]
tim=cf.jd-16000.  & data=cf.mnvel  & err=cf.errvel 
;pergram,tim,data,lowper=1.1, pmax=100.
;plots,[3.24, 3.24],[0,10],linesty=2,col=222,thick=4

wait,2
print,'Plotting the time series cleaned data' 
coef=poly_fit(tim,data,2)
plot,tim,data,ps=8,title='!6 Cleaned aCen B velocities', yra=[-70,50], /ysty
oplot,tim,poly(tim,coef),linesty=2,col=222,thick=4
oploterr, tim, poly(tim,coef), err

;wait,2
;print,'Plotting the data with the trend removed' 
;resid_vel=data-poly(tim,coef)
;plot,tim,resid_vel, ps=8, title='!6 Cleaned and detrended aCen B velocities', yra=[-20,20],/ysty
;oploterr, tim, resid_vel, err

;wait,2
;print,'Plotting pergram of detrended data'
;pergram,tim,resid_vel,lowper=1.03, pmax=100
;plots,[3.24, 3.24],[0,20],linesty=2,col=222,thick=4

save, sdev, f='sdev_'+tag+'.dat'
velplot,cf3[good],yra=[-50,50]
cf3=cf3[good]
cf1=cf3
save,cf1,cf3,f='vstbank/vst'+star+'_clean_'+tag+'.dat'


end ;pro 
