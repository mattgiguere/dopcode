function expand_wc, wc, wavout, nord=nord, wav=wav, npix=npix
nord = n_elements(wc[0,*])
x = indgen(nord)
z = where(wc[0,*] gt 0, nz)
new = wc*0
wavcond = keyword_set(wav) and keyword_set(npix) 
if wavcond then begin
    wavout = fltarr(npix, nord)
    xpix = indgen(npix)
endif
for i = 0, 3 do begin
    a = polyfit(x[z], wc[i,z], 4)
    new[i,*] = poly(x, a)
endfor
if wavcond then for i = 0, nord-1 do wavout[*, i] = poly(xpix, new[*,i])

return, new
end
