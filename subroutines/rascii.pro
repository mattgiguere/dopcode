pro rascii,vname,ncol,fname,skip=skip
;skip is the number of lines at the top of the file that are to be skiped
    fud='?'  &  crow=-1
	openr,u,fname,/get_lun
    if n_elements(skip) eq 1 then for n=1,skip do readf,u,fud
	vec=fltarr(ncol)  &  vname=fltarr(ncol,long(1.E5/ncol))
	while (eof(u) eq 0) do begin
	   readf,u,vec  &  crow=crow+1  &  vname(*,crow)=vec
    end  ;while
	close,u
	vname=vname(*,0:crow)
return
end
