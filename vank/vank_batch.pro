pro vank_batch

readcol,'starname_kt.txt',nmk,f='a10'
for i=0,n_elements(nmk)-1 do begin
	ff=findfile('/mir1/files/vdf*'+strcompress(nmk(i))+'*',count=num)
	if num gt 2 then begin
	  vank,strcompress(nmk(i),/remove_all),'f','vdf',fit_key=2,mct=1500
	end
end

stop
readcol,'starname_lt.txt',nm,f='a10'
for i=0,n_elements(nm)-1 do begin
	ff=findfile('/mir1/files/vdf*'+strcompress(nm(i))+'*',count=num)
	if num gt 2 then begin
	  vank,strcompress(nm(i),/remove_all),'f','vdf',fit_key=1, mct=1800
	end
end


end
