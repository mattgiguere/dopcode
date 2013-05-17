pro ck_vel_ha

restore,'vstbank/vst10700.dat'
x=where(cf3.errvel lt 1.25)  & cf3=cf3[x]

num=n_elements(cf3)
ha=fltarr(num)
vel=fltarr(num)
err=fltarr(num)

ind1=[0]
for i=0,num-1 do begin
   date=strmid(cf3[i].obnm,4,6)
   rdsk,hd,'/tous/mir7/iodspec/'+date+'/'+cf3[i].obnm,2
   ha_str=sxpar(hd,'HA',count=count)
   if count eq 1 and ha_str ne 'hour_angle' then begin
   	  ind1=[ind1,i]
      res=strsplit(ha_str,':')
      deg=strmid(ha_str,0,res[1]-1)
      min=strmid(ha_str,res[1],2)
      sec=strmid(ha_str,res[2],4)
      ha[i]=ten(deg,min,sec)
      if deg eq '-00' then ha[i]=ha[i]*(-1.0)
      vel[i]=cf3[i].mnvel
      err[i]=cf3[i].errvel
   endif
endfor
nn=n_elements(ind1)
ind=ind1[1:nn-1]

stop
x=where(ha ne 0.0)
ha=ha[ind]
vel=vel[ind]
err=err[ind]
i=sort(ha)  
ha=ha[i]  & vel=vel[i]  & err=err[i]
plot,ha,vel,ps=8,xtitl='!6 HA',ytitl='!6 Vel [m s!u-1!n]'

coef=poly_fit(ha,vel,3)
oplot,ha,poly(ha,coef),ps=8,col=90

stop

plot,ha,err,ps=8,xtitl='!6 HA',ytitl='!6 Error [m s!u-1!n]'
coef=poly_fit(ha,err,3)
oplot,ha,poly(ha,coef),ps=8,col=110
stop

end 
