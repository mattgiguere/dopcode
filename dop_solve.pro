;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_SOLVE, CHUNK, DOPENV, PARINFO, WT, PAR, SIGPAR, CHI,
; IODFLAG=IODFLAG, PLOT=PLOT
;
; substitute for mpfitfun in dop_marq
; this is the old Marcy/Butler/Valenti LM fitting routine
;
; INPUT:
; chunk
; dopenv
; parinfo
; wt
; pass
; iod_flag
; 
; OUTPUT:
; par
; sigpar
; chi
; rms 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION DOP_SOLVE,chunk,dopenv,parinfo,wt,pass,iod_flag=iod_flag, $
   fit=yfit, resid=resid, sigpar=sigpar,chi=chi,plot=plot

  npar = n_elements(parinfo.value)
  obchunk = (*chunk.sobs) 
  wiod = (*chunk.wiod)
  siod = (*chunk.siod) 
  wstar = (*chunk.wis)
  sstar = (*chunk.sis)
  osamp = dopenv.osamp
  numpix = dopenv.n_pix

 ; DEFINE PAR'S FROM THE PARINFO STRUCTURE
 ; par = [w0, disp, z, cont, ampl0, ampl(1-16)]
  par=fltarr(npar)   ; free parameters
  for i=0,npar-1 do par[i]=parinfo[i].value
  dstep=fltarr(npar) 
  for i=0,npar-1 do dstep[i]=parinfo[i].step
  w0 = par[0]
  dispobs = par[1]
  fltpar = where(parinfo.fixed ne 1,nfreepar)
  newpar = par 
  if ~keyword_set(iodflag) then iodflag=0
  spec = obchunk 

;Initialize some parameters.
  degf = n_elements(where(wt gt 0)) - npar - 2  ;degrees of freedom in fit
  c = 2.99792458d8                              ;also !c
  lambda = 0.3d0           		        ;init NLLS-gradient weight
  lamfac = 10.d0                                ;factor to change lambda
  idiag = indgen(nfreepar) * (nfreepar+1)	;indicies of main diagonal
  uvec = fltarr(nfreepar) + 1			;make unit vector

;Initial Synthetic Spectrum: yfit
  x=findgen(dopenv.n_pix)
  functargs = {chunk:chunk, dopenv:dopenv, pass:pass, iod_flag:iodflag}
  yfit = dop_fit(x, par, _EXTRA=FUNCTARGS)
  if n_elements(yfit lt 2) then stop

;Initial CHI-SQ
  resid = double(spec - yfit)                   ;compute reduced chi-sq
  wtres = wt * resid^2.0                        ;weighted residuals^2
  oldchisq = total(wtres) / degf	        ;chi-squared
  chisq=oldchisq


;Cross Correlate to refine input wavelength zero pt.
  index = where(fltpar eq 0,nind)               ;Is par(11) floating?
  IF nind gt 0 then begin                       ;Yes.  Refine it...
       xcorlb,spec,yfit,5,shft                  ;shft is pixel shift
       par[0] = par[0] - shft*dispobs           ;apply shft to wav zero pt
  ENDIF

;Initialize parameters for iteration
;;;  tpar = transpar(par,1,numpix)        ;New Transformed Params
  tpar = par  ; bypass tpar
  niter = 0
  iter = 0 
  chicount = 0
  maxiter = 100
  crit = 0.2                           ;Convergence Crit = 0.2 dstep = 2m/s)
  qq = where(dstep gt 0.)

  REPEAT BEGIN                                  ;ITERATION LOOP
     iter = iter + 1
     yfit = dop_fit(x, par, _EXTRA=FUNCTARGS)
     pder = dblarr(numpix, nfreepar)
     for n=0,nfreepar-1 do begin 
        ind=fltpar[n]
        dumpar = double(tpar)
        dumpar[ind] = tpar[ind] + dstep[ind]
;        par = transpar(dumpar, -1) 
        par = dumpar   ; no transpar
        hifit = dop_fit(x, par, _EXTRA=FUNCTARGS)
        if n_elements(hifit) lt 2 then stop
        pder[*,n] = (hifit - yfit) / dstep[ind]
        izero = where(pder[*,n] eq 0., nzero)
        if nzero ge 1 then pder[izero,n] = 1.d-50*median(yfit)/dstep[ind] 
     endfor
     if n_elements(yfit) le 2 then stop
     resid = double(spec-yfit)                  ;residuals
     wtres = wt * resid^2.0 
     oldchisq = total(wtres) / degf	        ;calculate chi-squared
     beta = double(resid*wt # pder) 
     alpha = double(transpose(pder) # (wt#uvec * pder))
     norm = sqrt(alpha(idiag) # alpha(idiag))  ;norm of diagonal elements
     array0 = alpha / norm		       ;normalized "alpha"
     REPEAT BEGIN                              ;LAMBDA LOOP
        chicount = chicount + 1 
        array = array0
        array(idiag) = 1.0 + lambda		;set LS vs. gradient search
        array = invert(array)			;invert array
        dfree = float(array/norm # transpose(beta)) ;parameter adjustments
        newpar(fltpar) = tpar(fltpar) + dfree	;try new free param values
;        newpar(0) = median([newpar(0),0.2,1.6]) ;central gaussian sigma limits
        yfit = dop_fit(x, newpar, _EXTRA=FUNCTARGS)
        if n_elements(yfit) le 2 then stop
        resid = double(spec - yfit)             ;residual of model fit
        wtres = wt * resid^2                    ;weighted residual squared
        chisq = total(wtres) / degf             ;calculate chi-squared
        lambda = lambda * lamfac                 ;assume fit worse;big lam => small step
     ENDREP UNTIL (chisq le oldchisq or lambda gt 1.d6)  ;this lambda gave smaller chi
     lambda = lambda/lamfac^2                  ;prepare lambda for next iter
     if lambda lt 1.d-4 then lambda = lambda*lamfac  ;avoid lambda too low
;    CONVERGED?
     delta = newpar - tpar                  ;change of params, this iter 
     if chisq le oldchisq then tpar = newpar 
     test= delta[qq]/dstep[qq]
     IF iter ge maxiter   or $              ;Convergence taking too long
        lambda gt 1.d6  then begin
           chi = 100.
           stop
     ENDIF
  ENDREP UNTIL total(abs(test)) lt crit ;end iter loop 0,100
  par = tpar
  chi = sqrt(chisq) 
  niter = fix(iter)

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


; BOMB:    ;emergency escape from a BOMB condition

RETURN, par 
END             ;end whole program
