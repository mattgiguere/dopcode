pro ck_vd, tag=tag


restore,'/mir7/bary/qbcvel.dat'
x=where(bcat.objnm eq '128620' and bcat.jd gt 14830. and bcat.obtype eq 'o',nx)
if ~keyword_set(tag) then tag = 'b'
openw,1,'missing_vd.txt'

for i=0,nx-1 do begin 
   vdnm=strcompress('vd'+tag+'128620'+bcat[i].obsnm,/rem)
   ff=file_search('/mir7/files_df/'+vdnm,count=count)
   if count eq 0 then printf,1,bcat[i].obsnm
endfor 

close,1

end
