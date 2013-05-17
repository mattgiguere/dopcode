pro ck_chi

ff=file_search('/tous/mir7/files_df/vdadgHR*',count=num)

for i=0,num-1 do begin
   restore,ff[i]
   if min(vd.fit eq 0.0) then print,'bad: ',ff[i]
   if i eq 0 then plothist,vd.fit,bin=0.1,xra=[0,10]
   if i gt 0 then plothist,vd.fit,bin=0.1, /overplot
   xx=where(vd.fit gt 4.0,nxx)
   if nxx gt 0 then print,vd[xx].ordt,vd[xx].pixt
endfor 

stop
end
