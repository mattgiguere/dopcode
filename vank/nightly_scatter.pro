pro nightly_scatter, tag=tag, archive=archive

;dfd='/tous/mir7/files_df/'
;if keyword_set(archive) then dfd=dfd+archive
;vank,'10700',tag,'vd'+tag,cf1,cf3,/ctio, dfd=dfd
vstpath='/Users/debrafischer/dop2/vank/vstbank/'
restore,vstpath+'vst10700.dat'
x=where(cf3.jd lt 16140.)  & cf3=cf3[x]
fjd=intarr(n_elements(cf3))

for i=0,n_elements(fjd)-1 do fjd[i]=fix(cf3[i].jd)
g=uniq(fjd)

fjd=fjd[g] 
fobnm=strarr(n_elements(fjd))
nrms=fltarr(n_elements(fjd))

for i=0,n_elements(fjd)-1 do begin
	x=where(fix(cf3.jd) eq fjd[i],nx)
	if nx ge 2 then begin
		nrms[i] = stddev(cf3[x].mnvel)
		fobnm[i] = cf3[x[0]].obnm
	endif
endfor 
	
yy=where(nrms eq 0.0,nyy) 
if nyy gt 0 then stop, 'drop some nights' 

xlo=where(nrms le 2.0,nxlo)
xbig2=where(nrms gt 2.0 and nrms le 3.0,nxbig2)
xbig3=where(nrms gt 3.0,nxbig3)
print,'the number of nights where the rms scatter is le 2 m/s: ', nxlo
if nxbig2 gt 0 then print, 'the nights where rms scatter is 2.0-3 m/s: '
for j=0,nxbig2-1 do print, j+1,' ',fjd[xbig2[j]],' ',fobnm[xbig2[j]],' ',nrms[xbig2[j]]
if nxbig3 gt 0 then print, 'the nights where rms scatter is greater than 3.0 m/s: '
for j=0,nxbig3-1 do print, j+1,' ',fjd[xbig3[j]],' ',fobnm[xbig3[j]],' ',nrms[xbig3[j]]

stop
end ;pro
