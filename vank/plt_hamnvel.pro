pro plt_hamnvel 

restore,'/home/fischer/dop2/vank/vstbank/vst161797.dat'
iodpath='/mir1/iodspec/'

num=n_elements(cf3.obnm)
ha=fltarr(num) 

for i=0,num-1 do begin
   rdsk,hd,iodpath+cf3[i].obnm,2
   ha_tmp=strcompress(sxpar(hd,'HA'),/rem)
;   a1=strpos(ha_tmp,':')
   rah=strmid(ha_tmp,0,2)
;   ha_tmp2=strmid(ha_tmp,a1+1,2)
;   a2=strpos(ha_tmp2,':')
   ram=strmid(ha_tmp,3,2)
   ras=strmid(ha_tmp,6,4)

   ha[i]=ten([rah,ram,ras])

print,ha[i]
print,''
end

plot,ha,cf3.mnvel,ps=8


stop

end
