function dop_loc, x, par, _EXTRA=FUNCTARGS

;fischer oct 2009
; x is the findgen(80) array of pixels 
; par contains the parameters needed to stretch and 
;    shift the FTS atlas onto the observed chunk.
; par = [w0, disp, shift] 

c_light=2.99792458d8

extra=functargs
wiod=extra.wiod   ;I2 wavelengths
siod=extra.siod   ;I2 spectrum
sobs=extra.btmp_ch ; observed 80-pixel I2 spectrum
wid=extra.wid

ip=psf(wid)
synth_iod=convol(siod,ip,/edge_truncate, /normalize)

wobs = par[0] + par[1]*(x+par[2])
rebin, wiod, synth_iod, wobs, synth_new

scale = sobs/synth_new
contf,scale,c_scale,sbin=10, nord=1, frac=0.3
syn_fit = par[3]*c_scale*synth_new

return, syn_fit

end

