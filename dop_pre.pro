;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_PRE
;
; CALLED FROM: DOP_MARQ
;
; CALLS TO: DOP_FIT 
; 
; PURPOSE:  
;   Drives Levenberg-Marquardt search to refine wavelength soln
;   and Doppler z.  Central Gaussian width is only other free parameter
;
; PROCEDURE:
;   1) weight the pixels 
;   2) set up the parinfo structure for free parameters
;   3) call dop_pre_fit to create a synthetic spectrum
;   4) calculate a reduced chi-sq fit
;   5) return the best model free parameters to dop_chunk_setup
; 
; INPUTS: 
;   CHUNK
;      sobs: observed spectrum
;      siod: FTS spectrum;        wiod: FTS wavelengths
;      sis: intrinsic spectrum;  wis: wavelengths of intrinsic
;      psfpix, psfsig
;      chunk.free_pars
;   DOPENV
;      osamp      
; 
; OUTPUTS:
;   MODEL
;      chunk.free_pars (updated) 
;   
; OUTSTANDING ISSUES: 
; the psf peak is not centered on the 0 pixel b/c the 
; side gaussians kick in.  this seems non-physical. 
;
; Written by: Debra Fischer, SFSU, Dec 2007
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_PRE, chunk, npre_free, dopenv, demo=demo, lm=lm, pfilt=pfilt, avg=avg

    c_light = dopenv.c_light
    sobs = (*chunk.sobs)

    if ~keyword_set(demo) then demo=0

  ; WEIGHT THE PIXELS
    noise_flat = 0.005   ;typical s/n in flatfield
    wt = (1./sobs) / (1. + sobs*noise_flat^2)
    xbad=where(wt lt 0. or wt gt 5.*median(wt),nbad)   
    if nbad gt 0 then wt[xbad]=0.0
    wt = wt*pfilt   ; zero out the bad pixels 
	sig = sqrt(*chunk.sobs) ;only used in plotting

  ; INFORMATION AND LIMITS FOR FREE PARAMETERS: PARINFO STRUCTURE
  ; THIS IS EXPLICITLY THE PRE-FIT, SO ALWAYS PASS=0
    nfree=npre_free
    struct = (*chunk.free_par)[0]  ;dop_bearings only fills in pass[0]

  ; SETUP THE PARINFO STRUCTURE
    parinfo={value: 0.0d,       $ ; double precision  
             fixed: 0,          $
             limited: [0,0],    $    ; use with caution 
             limits: fltarr(2), $
             parname: '?',      $
             step: 0.0d,        $
             relstep: 0.00,     $
             mpside: 0}              ; 0, 1, -1, 2   
    parinfo=replicate(parinfo, npre_free)

    parinfo[0].value = double(struct.wcof[0])    ;wavelength
    parinfo[0].step = 0.001
    
    parinfo[1].value = double(struct.wcof[1])    ;dispersion
    parinfo[1].step=0.0001
;    if dopenv.n_wcof eq 3 then begin
;    	ij=1
;    	parinfo[1+ij].value=double(struct.wcof[2])
;     	parinfo[1+ij].step = 0.0001
;    	parinfo[1+ij].limited = [1,1]
;    	parinfo[1+ij].limits = [0.96*parinfo[1].value, 1.04*parinfo[1].value]/10.
;    endif 
    if dopenv.n_wcof eq 2 then ij=0
   	
    parinfo[2+ij].value = double(struct.z*c_light)  ;Doppler vel (including BC)
    parinfo[3+ij].value = struct.offset     ;continuum scaling
;jan 7,2013 fischer 

    ;print,parinfo[4+ij].value
     ;  stop
  ; PAR[2]: DOPPLER SHIFT        
    parinfo[2+ij].step = 5.   ;1.e-6
    if dopenv.iss_nm eq 'iod' then parinfo[2+ij].fixed=1

  ; PAR[3]: SCALE CONTINUUM (THIS USED TO BE SCATTERED LIGHT) 
    parinfo[3+ij].step = 0.01
    ;parinfo[3+ij].limited = [1,1]
    ;parinfo[3+ij].limits = [0.95, 1.05]

;jan 7,2013 fischer 
  ; PAR[4]: CENTRAL GAUSSIAN WIDTH
  if nfree eq 5 then begin
     parinfo[4+ij].value = struct.amp[0]   ;central Gaussian width - free par
     parinfo[4+ij].step=0.001 
     parinfo[4+ij].limited = [1,1]
     parinfo[4+ij].limits = [struct.amp[0]*(1.0-0.2),struct.amp[0]*(1.0+0.2)]
;print,'dop_pre psf limits: ',parinfo[4+ij].limits
  endif
;print,'limiting the width of the central gaussian: ',parinfo[4+ij].limits
    dof=n_elements(sobs) - nbad - npre_free ; weight=0 or fixed
    
  ; THINGS CAN ONLY GET BETTER FROM HERE... 
    pass=0
    functargs = {chunk: chunk, dopenv: dopenv, pass:pass}
    x = findgen(dopenv.n_pix)    
    y = sobs

  ; NLLS FITTING - TURN THE (TEMPLATE*IODINE) INTO THE OBSERVATION
         newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, avg=avg, $
                         functargs=functargs, maxiter=300,/nan, $
                         yfit=syn_fit, weight=wt, status=status, /quiet)

if n_elements(syn_fit) eq 0 then stop
  ;quick cross correlation
;	xcorlb,sobs, syn_fit, 10, shft
;	parinfo[0].value = parinfo[0].value - shft*parinfo[1].value
;	newpar[0] = newpar[0] - shft*newpar[1]
	
;	plot,sobs,ps=8
;	oplot,syn_fit,col=222
;	print,shft
;	stop
	
  ; error checking (rebin wavelength out of range)  17nov2011 dfischer 
       if n_elements(newpar) eq n_elements(parinfo.value) then parinfo.value = newpar 
       if n_elements(newpar) eq 1 then begin 
          tries = 0L
          wavperturb = 0.005d
          par0 = parinfo[0].value
          repeat begin ;added 20120126
          parinfo[0].value= par0 + randomn(seed)*wavperturb
          newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, avg=avg, $
                          functargs=functargs, maxiter=200, $
                          errmsg=errmsg, /iterstop, $
                          bestnorm=bestnorm, perror=perror, $
                          yfit=syn_fit, weight=wt, status=status,/quiet)
          tries++
          endrep until (n_elements(newpar) eq n_elements(parinfo.value) OR tries eq 1d3)
		if n_elements(newpar) lt n_elements(parinfo.value) then newpar=parinfo.value
		parinfo[0].value=parinfo[0].value+randomn(seed)*0.001
        if n_elements(newpar) eq n_elements(parinfo.value) then parinfo.value = newpar
      endif
      
   ; UPDATE THE CHUNK INITIAL PARAMETERS FROM THIS PRE-PASS 
      (*chunk.smod)=syn_fit
   ; TOUCH UP THE WAVELENGTH SOLUTION FOR PASS 1 
   ; SMOD IS FROM DOP_PRE 
;      xcorlb,(*chunk.sobs),(*chunk.smod),15,shift  ;2013jan12 comment out
;      newpar[0] = newpar[0] - shift*newpar[1]     ;2013jan12 comment out
      (*chunk.free_par)[0].wcof[0] = double(newpar[0])
      (*chunk.free_par)[0].wcof[1] = double(newpar[1])
      if ij eq 1 then (*chunk.free_par)[0].wcof[2] = double(newpar[2])
      (*chunk.free_par)[0].z = double(newpar[2+ij]/c_light)
      (*chunk.free_par)[0].offset = newpar[3+ij]
 ;jan 7,2013 fischer 
      if nfree eq 5 then (*chunk.free_par)[0].amp[0] = newpar[4+ij]  ; really the gaussian width, not ampl

      red_chi_sq = total(  wt * (syn_fit - sobs)^2 )
      red_chi_sq = red_chi_sq / dof
      chi = sqrt(red_chi_sq) 
      
      chunk.fit[0] = chi
      ;print,'Pass 0: Red chi ',chi, ' FWHM: ',newpar[4]
      if demo eq 1 then begin 
      ; COLORS FOR PLOTS
         loadct,39, /silent
         !p.background=255
         !p.color=1

;         ;print,'Pre-pass reduced chi-sq: ', chi
;
       ; FIRST PANEL: PSF
         !p.charsize=1.8
         !x.charsize=1.6
         !y.charsize=1.6
         !p.multi=[0,1,3]
         !x.omargin=[6,2]
         !y.omargin=[2,2]
         ip=dop_pre_psf(newpar[4],dopenv=dopenv)
         plot, ip, col=1, title='IP', /xsty,xra=[40,80], /ysty, yra=[0,1]
         plots, [60,60], [0,1]

       ; SECOND PANEL: DECONV SPECTRUM AND FTS IODINE
         plot, (*chunk.wis), (*chunk.sis),col=1,title='!6 ISS and I2', $
               /xsty, yra=[0,1.1],/ysty, yticks=2,$
               ytickv=[0, 0.5, 1.0], ytickname=['0','0.5','1.0']
         oplot, (*chunk.wiod), (*chunk.siod), col=100

         wav_coef = (*chunk.free_par)[0].wcof
         wobs = poly(findgen(dopenv.n_pix), wav_coef)

       ; THIRD PANEL: MODEL AND OBSERVED SPECTRUM
         plot, wobs, (*chunk.sobs), xra=[min(wobs),max(wobs)], /xsty, yticks=3,$
               /ysty, xtit='!6 Wavelength', $
               titl='!6 Pre Model (blue) and Obs Spectrum (black)'
         !p.color=70
         oploterr, wobs, syn_fit, sig
         oplot, wobs, syn_fit, col=70
         !p.color = 1
      endif

   refine_ch = chunk
   return, refine_ch

end

