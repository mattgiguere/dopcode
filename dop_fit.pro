;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_FIT
;
; initial coding: Debra Fischer, Dec 2007
;
; CALLED FROM: DOP_MARQ
;
; PURPOSE:  
;   Creates the synthetic observation with free parameters 
;   (PSF description, pixel wavelengths, doppler shift)
;   driven by a Levenberg-Marquardt algorithm
;
; PROCEDURE:
;   1) multiply FTS and DST
;   2) convolve with PSF
;      - free parameters: PSF description (PASS 1)
;      - Doppler shift, 
;        wavelength solution, continuum offset
; 
; INPUTS: 
;   X: independent variable (pixels)
;   PARINFO: initial values 
;
; -EXTRA:
;   CHUNK: observed spectrum (with iodine)
;   DOPENV: structure with paths, psf descrip etc
; 
; OUTPUTS:
;   PAR: best fit free parameters for mpfit
;
; OUTSTANDING ISSUES: 
; put in a check - if the velocity changes by too much too little, patch
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_FIT, X, PAR, _EXTRA=FUNCTARGS

COMMON MPFIT_ERROR, ERROR_CODE 

  ; UNPACK THE VARIABLES 
    extra=functargs
    dopenv = extra.dopenv 
    osamp = dopenv.osamp
    avg = dopenv.avg
    chunk = extra.chunk
    pass = extra.pass
    c_light = dopenv.c_light
    
   ; print, par
   ; print, 'n_elements in par in dop_fit: ', n_elements(par)

  ; WAVELENGTH OF THE OBSERVATION
    disp = par[1]
    coef = [par[0], par[1]]  
    dshift = par[2]/c_light
    wobs = poly(x, coef)     ; wavelength of observed spectrum (from free pars)

  ; DOPPLER SHIFT THE TEMPLATE 
    wis_shft = double((*chunk.wis)) + double((*chunk.wis))*dshift 
    ;lambda_new = (1+ v/c)*lambda

  ; FINE WAVELENGTH SCALE FOR THE MODEL
    w0fine = max( [(*chunk.wiod)[0], wis_shft[0]] ) + 0.5*disp
    w1fine = min( [max((*chunk.wiod)), max(wis_shft)] ) - 0.5*disp
    disp_fine = double(par[1] / osamp)
    nfine = fix( (w1fine - w0fine) / disp_fine )  ;extra padding
    wav_fine = w0fine + dindgen(nfine)*disp_fine  

;!p.multi=[0,1,3]
;plot,(*chunk.wiod),(*chunk.siod),/xsty
;plot,wis_shft,(*chunk.sis),/xsty 
;plot,wobs,(*chunk.sobs),/xsty
;!p.multi=[0,1,1]

error_code=0
  ; WAV_FINE RANGE SHOULD BE GREATER THAN THAT OF 
  ; THE OBSERVATION, BUT A SUBSET OF THE IODINE AND TEMPLATE
    if wav_fine[0] gt wobs[0] then error_code=-1              ; doesn't cover observation
    if max(wav_fine) lt max(wobs) then error_code=-2          ; doesn't cover observation
    if wav_fine[0] lt (*chunk.wiod)[0] then error_code=-3     ; outside the iodine range
    if max(wav_fine) gt max((*chunk.wiod)) then error_code=-4 ; outside the iodine range
    if wav_fine[0] lt wis_shft[0] then error_code=-5          ; outside the template range
    if max(wav_fine) gt max(wis_shft) then error_code=-6      ; outside the template range
if error_code ne 0 then stop

  ; SPLINE IODINE AND DSST ONTO OVERSAMPLED DSST WAVELENGTH SCALE
    siod_fine=dspline((*chunk.wiod),(*chunk.siod),wav_fine)
    sis_fine=dspline(wis_shft,(*chunk.sis),wav_fine)

  ; DECONVOLVED MODEL SPECTRUM
    dss = siod_fine*sis_fine
    
  ; CONVOLVE WITH IP, BIN TO OBSERVED WAVELENGTH SCALE,
  ; REGISTER THE OBSERVED AND SYNTHETIC SPECTRA ON THE Y-AXIS
;jan7, 2013 fischer
    if pass eq 0 and avg eq 0 then ip = dop_pre_psf(par[4], dopenv=dopenv)
	if pass eq 0 and avg eq 1 then ip = chunk.ip[*,1]  ; 4 free pars - not psf
    if pass eq 1 then begin
    	if dopenv.psfmod eq 'gaussian' then ip = dop_psf(par[4:n_elements(par)-1],dopenv=dopenv)
    	if dopenv.psfmod eq 'bspline' then ip = dop_psf_bspline(par[5:n_elements(par)-1], dopenv=dopenv, cntr=dopenv.psfcntr)
    endif   ; pass eq 1
    if pass eq 2 then ip = dop_psf_smooth(chunk, chunk.ordob, chunk.pixob, $
                                          dopenv=dopenv, chunk_arr=extra.chunk_arr)
    synth_spec=convol(dss, ip, /edge_truncate, /normalize)
    rebin, wav_fine, synth_spec, wobs, synth_new

  ; MATCH FLUX FOR THE OBSERVED AND SYNTHETIC SPECTRA 
    scale = (*chunk.sobs) / synth_new
    contf,scale,c_scale, sbin=20, nord=3, frac=0.5
    syn_fit = par[3]*c_scale*synth_new
    
    return, syn_fit

end
