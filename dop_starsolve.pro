Function STARSOLVE,par,sigpar,chi,plot=plot

common spectra,obchunk,wiod,siod,wstar,sstar,w0,inpr,osamp,wt,$ 
               dstep,npar,fltpar,keck
common psfstuff,psfsig,psfpix,obpix

;PURPOSE:
;  Find free parameters, PAR, to fit synthetic to Observed Spectrum.
;
;INPUT:
;  PAR     fltarr(15)   Input Guesses of Free Parameters.
;    PAR (0:10) = IP parameters (see gpjv.pro or gpfunc.pro)
;    PAR(11) = wavelength zero point
;    PAR(12) = cZ
;    PAR(13) = dLambda/dpix
;  OBCHUNK fltarr(40)  Observed spectrum --- #photons in each pixel
;  WIOD    fltarr      Wavelengths for fts spectrum (iodine)
;  SIOD    fltarr      Fts spectrum (iodine) --- finely sampled spectrum
;  WSTAR   fltarr      Wavelengths of stellar spectrum (deconvolved)
;  SSTAR   fltarr      Stellar Spectrum (deconvolved)
;  W0      int         integer portion of wavelength of chunk
;  OSAMP   int         # of sub-pixels per original pixel (=4)
;  WT      fltarr      Pixel Weights in OBCHUNK:  1/eps^2 from photon stat
;  DSTEP   fltarr(15)  Differential increments in free parameters (~10m/s)
;  NPAR    int         # of parameters
;  FLTPAR   fltarr      Indices of 15 param's that are free to float.
;  KECK                Flags use of Keck/HIRES spectrum
;  TRACE:  optional output - higher values of trace ==> more screen diagnositcs
;
;OUTPUT:
;   PAR    fltarr(15)   Values of free parameters that best fit OBCHUNK
;   SIGPAR: Uncertainties in IP parameters
;   CHI:    chi-sq for fit
;   RMS:    optional output - rms to fit
;----------------------------------------------------------------
;
;Check for absurd input parameters
  if par(11) lt -1. or par(11) gt 2. or par(13) le 0. then begin
      print,'STARSOLVE: Abort because of absurd input wavelength scale.'
      print,'par(11)=',par(11)
      print,'par(13)=',par(13)
      chi = 100. 
      return,-1
  end
;
;Initialize some parameters.
  dispobs = par(13)                             ;dispersion - dlambda/dpixel
  numpix = n_elements(obchunk)                  ;# of pixels in OBCHUNK (40)
  spec = obchunk                                ;rename "obchunk": "spectrum"
  newpar = par                                  ;array for new par's
  degf = n_elements(where(wt gt 0)) - npar - 2  ;degrees of freedom in fit
  c = 2.99792458d8                              ;also !c
  lambda = 0.3d0           		        ;init NLLS-gradient weight
  lamfac = 10.d0                                ;factor to change lambda
  idiag = indgen(npar) * (npar+1)		;indicies of main diagonal
  uvec = fltarr(npar) + 1			;make unit vector

;Initial Synthetic Spectrum: yfit
  tpar = transpar(par,1)                        ;Transform params, dsteps
  yfit = starham(tpar)	                        ;synthetic "fit" spectrum

;Check for serious flaws in synthetic spectrum or weights
       IF (n_elements(yfit) le 1) or (degf le 1) or $
	(n_elements(where(yfit gt 0)) lt 2) then begin
          print,'STARSOLVE BOMB:' 
          print,'  n_elements(yfit)=',n_elements(yfit)
          print,'  degf=',degf
          print,' '
          chi=100. & goto,BOMB
       ENDIF
;
;Initial CHI-SQ
  resid = double(spec - yfit)                   ;compute reduced chi-sq
  wtres = wt * resid^2.0                        ;weighted residuals^2
  oldchisq = total(wtres) / degf	        ;chi-squared
  chisq=oldchisq
;

;Cross Correlate to refine input wavelength zero pt.
  index = where(fltpar eq 11)                    ;Is par(11) floating?
  IF index(0) ge 0 then begin                   ;Yes.  Refine it...
       xcorlb,spec,yfit,5,shft                  ;shft is pixel shift
       par(11) = par(11) - shft*dispobs         ;apply shft to wav zero pt
  ENDIF

;Initialize parameters for iteration
  tpar = transpar(par,1)                        ;New Transformed Params
  iter = 0 
  maxiter = 100
  crit = 0.2                    ;Convergence Crit = 0.2 dstep = 2m/s)

  REPEAT BEGIN                                  ;ITERATION LOOP
     iter = iter + 1
;        if par(0) lt 0. then par(0) = 0.
;        if par(2) lt 0. then par(2) = 0.
;     ineg = where(tpar(0:4) lt 0., nneg)       ;Demand positive IP pars
;     if nneg ge 1 then tpar(ineg) = 0.          ;pin at zero
     yfit = starham(tpar,pder)                  ;First synthetic spectrum
     if n_elements(yfit) le 2 then begin
	print,'STARSOLVE: BOMB in starham' & chi=100. & GOTO,BOMB
     endif
     resid = double(spec-yfit)                  ;residuals
     wtres = wt * resid^2.0 
     oldchisq = total(wtres) / degf	        ;calculate chi-squared
     beta = double(resid*wt # pder) 
     alpha = double(transpose(pder) # (wt#uvec * pder))
     norm = sqrt(alpha(idiag) # alpha(idiag))  ;norm of diagonal elements
     array0 = alpha / norm		       ;normalized "alpha"
     REPEAT BEGIN                              ;LAMBDA LOOP
        array = array0
        array(idiag) = 1.0 + lambda		;set LS vs. gradient search
        array = invert(array)			;invert array
        dfree = float(array/norm # transpose(beta)) ;parameter adjustments
        newpar(fltpar) = tpar(fltpar) + dfree	;try new free param values
        newpar(0) = median([newpar(0),0.2,1.6]) ;central gaussian sigma limits
;        if par(0) lt 0. then par(0) = 0.
;        if par(2) lt 0. then par(2) = 0.
;        ineg = where(newpar(0:4) lt 0., nneg)   ;Demand positive IP pars
;        if nneg ge 1 then newpar(ineg) = 0.     ;pin at zero
;       par = transpar(newpar,-1)
;       forma = '(I6,8F10.5)'
;       print,format=forma,iter,tpar(0),tpar(3),tpar(4),tpar(5),tpar(6),tpar(7),tpar(8),sqrt(chisq)
;       forma = '(I6,9F10.5)'
;       print,format=forma,iter,dfree(0),tpar(0),tpar(3),tpar(4),tpar(5),tpar(6),tpar(7),tpar(8),sqrt(chisq)
;       print,'lambda=',lambda
        yfit = starham(newpar)                  ;compute synthetic profile
        if n_elements(yfit) le 2 then begin
           print,'STARSOLVE: BOMB in starham' & chi=100. & GOTO,BOMB
        endif
        resid = double(spec - yfit)             ;residual of model fit
        wtres = wt * resid^2                    ;weighted residual squared
        chisq = total(wtres) / degf             ;calculate chi-squared
        forma = '(a6,6F10.5)'
        lambda = lambda * lamfac                 ;assume fit worse;big lam => small step
     ENDREP UNTIL (chisq le oldchisq or lambda gt 1.d6)  ;this lambda gave smaller chi
     lambda = lambda/lamfac^2                  ;prepare lambda for next iter
     if lambda lt 1.d-4 then lambda = lambda*lamfac  ;avoid lambda too low
;    CONVERGED?
     delta = newpar - tpar                  ;change of params, this iter 
     tpar = newpar                          ;save these current par's
     dvel = delta(11:13)/dstep(11:13)       ;dpar in units of dstep (10 m/s) velocity
     IF iter ge maxiter   or $              ;Convergence taking too long
        lambda gt 1.d6  then begin
           chi = 100.
           goto,BOMB
     ENDIF
  ENDREP UNTIL total(abs(dvel)) lt crit          ;end iter loop 0,100

;print out some important parameters
;  pat = transpar(newpar,-1)
;  forma = '(I6,8F10.5)'
;  print,format=forma,iter,pat(0),pat(1),pat(2),pat(3),pat(4),pat(11),pat(13),sqrt(chisq)
;
par = transpar(tpar,-1)     ;TRANSFORM PARAMS BACK
chi = sqrt(chisq)
;
;
;PLOTTING SECTION
  if n_elements(plot) eq 0 then plot = 0		;default: no plots
  if plot eq 1 then begin
   !p.charsize=1.
   !x.charsize=1
   !y.charsize=1
    xwav = dindgen(numpix)*par(13)+par(11)+w0          ;linear wavelength scale
    xtit='!6Wavelength (A)'
    ytit='!6Residual Intensity'
    tit='!6 .+* : Obs Spec     __ : Synth Spec'
    fac = max(spec)
    plot,xwav,spec/fac,psym=7,/yno,xtit=xtit,ytit=ytit,tit=tit,symsize=1.2
    oplot,xwav,yfit/fac,thick=1.8,co=151
    wait,1
    xip = (findgen(120)-60.)/4.
    if n_elements(inpr) lt 2 then plot,xip,gpfunc(xip,par),ps=8,xr=[-4,4] $
       else plot,xip,inpr,ps=8,xr=[-4,4]
  endif

BOMB:    ;emergency escape from a BOMB condition

RETURN,yfit
END             ;end whole program
