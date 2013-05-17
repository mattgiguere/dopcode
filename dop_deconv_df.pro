;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_DECONV_DF
; 
; CALLED FROM DOP_DSST
; PURPOSE:
; 	Creates dspec, a function that describes the deconvolved spectrum
;   driven by mpfitfun, a Levenberg-Marquardt algorithm
; 
; PROCEDURE:
;
;
; INPUTS: 
; 
;
; _EXTRA:
;
;
;
; OUTPUTS:
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_DECONV_DF, par, xfine=xfine, spec=osamp_spec, psf=psf, node=node, $
                        sig=sig, newspec=spec, ngrid=ngrid, resid=diff, trim=trim, $
                        mov=mov

  spec = osamp_spec
  wid = par[0]
  p = par[1:ngrid]
  offset = par[ngrid+1:*]
  npar = n_elements(p)
  npix = n_elements(xfine)

;  plot, spec

  for i = 0, npar-1 do begin
     dist = xfine - node[i]     ; xe - node position
     tau = -alog(spec)
     amp = dc_gauss(dist, [p[i], offset[i], wid] ) + 1.
     spec = exp(-tau*amp)       ; < 1
  endfor

  num_conv, spec, psf, newconv   ;spec is the deconvolved spectrum
  diff = newconv - osamp_spec
  dspec = spec
  
  if keyword_set(mov) then begin
  	o = 0.2
    plot, xfine, osamp_spec+o, yr=[-0.2, 2], /ys, ps=8, syms=0.5
    oplot, xfine, newconv+o, co=155
    oplot, xfine, diff*5, ps=8, syms=0.5
    oplot, xfine, osamp_spec, co=0
    oplot, xfine, dspec, col=222
    hline, 1.4, lines=2
    oplot, node+offset, node*0+1.4+p*0.1, ps=8, co=155
    plots, median(xfine)+[-1,1]*wid, [0.2, 0.2]
endif

resid=diff[trim:npix-trim-1]/sig

return, resid 

end
