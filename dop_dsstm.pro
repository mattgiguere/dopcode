pro dsst_m,star,vdiod,filter,cf,dsst,inp_vdiod=inp_vdiod, $
	   keck=keck,sigpsf=sigpsf,pixpsf=pixpsf,plot=plot, $
           label=label

; This routine computes deconvolved spectra from observed "templates".
; The result in put into a data structure: DSST
;
; The routine uses the Nat. Solar Observatory Solar Flux spectrum
;   as the first "guess" for the stellar deconvolved spectrum.  It
;   modifies the NSO spectrum in 5 ways:
;  1.  Doppler shifts the NSO
;  2.  Artificially deepens all lines in a chunk using enhanced Optical Depth.
;  3.  Artifically broadenins it (V sini)
;  4.  Artificially brooadens lines with a Lorentzian (gamma is free)
;  5.  Add Gaussians to the NSO to correct deficiencies leftover.

;
;STAR:   (input)  Template star 2D array, (pixel_length, num_orders)
;VDIOD:  (input)  Any appropriate template "vd" structure
;FILTER: (input)  Array specifying the flawed pixels.  See filt6.dat, etc.
;cf:     not used
;DSST:   (output) Deconvolved star structure
;inp_vdiod (keyword: on/off) if invoked this uses the "vdiod" as the
;	   vd to carry out the deconvolution, and ignores the cf.
;         
common template,wavtemp,temp,wavfine,sunfine,psf,dstep
common psfstuff,psfsig,psfpix,obpix
common gaustuff,gauwav,gauwid

;Initialize parameters
 par = dblarr(70)
 c = 2.9979246d8

; len is no longer hardwired, but is read from input VD (vdiod), 15Jun98 PB
; len = 40          ;n_elements(sp)
 splen = n_elements(star(*,0))        ;pixels per order
 osmp = 4         ;oversampling
 gauwid = 1.2      ;Width of added Gaussians in original pixels
 par(1) = 1.0      ;Initial Tau Factor  (1.0 ---> no change)
 par(2) = 1.0      ;Vsini  (km/s)
 par(3) = 0.5      ;Gamma for Lorentzian
 psfsig=sigpsf
 psfpix=pixpsf
 vd=vdiod  &  oldvd=vdiod
 print,' '

 ;Choose chunk
 ;chunk = 101
; read,'chunk:',chunk
; n=chunk

;Deal with both new and old style VD's
tagnam=tag_names(vdiod)                             ;VD tag_names
ordindex=first_el(where(tagnam eq 'ORDER'))         ;order_index
if ordindex eq -1 then begin
   ordindex = first_el(where(tagnam eq 'ORDT'))     ;ORDER or ORDT?
end

contlev = 9.e4                         ;arb. continuum level
loord = min(vd.(ordindex))             ;minimum order
hiord = max(vd.(ordindex))             ;maximum order
imnorm = max(star(*,loord:hiord))      ;template normalization constant
normstar = contlev * (star/imnorm)     ;normalize template

if n_elements(label) ne 1 then label='vd' else label=strtrim(label,2)
if n_elements(dfn) ne 1 then begin              ;"Files" directory
   dfn='/mir1/files/'                           ;Berkeley files disk
endif else dfn=strtrim(dfn,2)

d_pix=70  &  d_ord=3  &  orddst = 15  ;Pixel & order ranges for PSF averaging

if n_elements(keck) ne 1 then keck=0  ;KECK OBSERVATION?
if keck eq 1 then begin                     
   d_pix = 100  &  d_ord = 1       ;Pixel and order ranges for PSF averaging
   orddst = 65.
   dfn = '/mir3/files/'            ;Set "files" directory
   print,'Keck observation, Averaging: delta_pix = 100,  del_ord=1 '
endif                               ;End Keck Section

if keyword_set(inp_vdiod) then begin
   print,'Using input "VDIOD" for deconvolution.' 
   end else begin
      n_vd=n_elements(cf)
      for m=0,(n_vd-1) do begin
        vdfn=dfn+label+'_'+strtrim(cf(m).obnm,2)
        restore,vdfn
        if m eq 0 then dum=vd else dum=[dum,vd]
      endfor ;m
      vd=dum
   endelse

last = n_elements(vd)-1                         ;index of last vd row
vdset = first_el( where(vd.(ordindex) eq vd(last).(ordindex) $
	 and vd.pixt eq vd(last).pixt) )        ;1st occurence of last chunk

n_lines = vdset+1                               ;# of unique line sets
;  New code to deal with arbitrary chunk sizes (40 or 50 pixel), 15Jun98, PB
len=vd(0).npix                                  ;assume vd.npix are all the same
dstlen = 256                                    ;for 40 pixel chunks
if len gt 45 then dstlen=296                    ;for 50 pixel chunks
if len eq 40 then $  ;assume either 40 or 50 pixel chunks, 15Jun98, PB
  dum =  {ddst,ordt:0,pixt:0,w0:double(0.),w1:double(0.),weight:0.,dst:fltarr(dstlen)} $
    else dum =  {ds50,ordt:0,pixt:0,w0:double(0.),w1:double(0.),weight:0.,dst:fltarr(dstlen)}
dsst = replicate(dum,n_lines)                   ;Define dsst structure, all rows
dsst.pixt = vd(0:n_lines-1).pixt & dsst.ordt=vd(0:n_lines-1).(ordindex) ;set pixls,ords


xip =  findgen(30*osmp+1)/osmp-15.           ;x-scale for IP
ordr=-1

if n_elements(index) gt 0 then n_lines=n_elements(index)

if n_lines gt 1 then sm_wav,vd,wc,plot_key=plot else wc=0  ;smooth wavelengths
   
   fwhm = 50 & nord = 4        ;for 1st CONT when length > 1000 pixels
IF splen lt 1000 then fwhm = 25 & nord = 4   ;for CONT when length < 1000

FOR m=0,(n_lines-1) do begin               ;Cycle through chunks
   n=m
   if vd(n).(ordindex) ne ordr then print,' '
   ordr = vd(n).(ordindex)                        ;order of current chunk
   place = vd(n).pixt                        ;pixel of current chunk
   ind = where(vd.pixt gt (place-2) and vd.pixt lt (place+2) $ ;vd indices of
		     and vd.(ordindex) eq ordr,nind)           ;same pix, ord
   ;nind = number of concatinated vd's
   wd = vd(ind)                              ;subset of vd, at same pix,ord
   IF nind gt 2 then begin
     medfit = median(wd.fit)                   ;median of fit
     if max(wd.fit) gt (1.5*medfit) then $   ;toss high fit, if too high
       wd = wd(where(wd.fit lt max(wd.fit)))
   ENDIF
   wt = float(wd.npix)/(wd.fit^2.)          ;set up weigths (npix/fit^2)
   wt = wt/total(wt)                        ;normalize weights
   dispersion = vd(n).wcof(1)               ;weighted mean dispersion
   if n_elements(wc) gt 1 then begin
      lambda=first_el(double(poly_fat([place],wc(*,ordr)))) 
   end else begin
    ind=where(vd.(ordindex) eq vd(n).(ordindex) and vd.pixt eq vd(n).pixt)
    lambda = double(vd(n).w0)+total(vd(ind).wcof(0)*vd(ind).fit)/total(vd(ind).fit)
    dispersion = double(total(vd(ind).wcof(1)*vd(ind).fit)/total(vd(ind).fit))
   endelse

   scat = total( wd.scat*wt )               ;weighted mean scattered Lt
   lo = max([0,place-30])                   ;lo pxl  (pixel - 30)
   slop = 100                               ;Total Length padded to 100 pxl
   hi = min([lo+slop-1,splen-1])            ;hi pxl (middle has 40 pxls)
   if hi eq (splen-1) then lo=splen-slop    ;treat case of chunk at end
   filt=reform(filter(lo:hi,ordr))          ;filter
   prop_filt,filt,/zero                   
   filt(0:12) = 0. &  filt(slop-13:slop-1)=0. ;Give no weight to ends
   xxnow = findgen(slop) + (lo-place)       ;[pix-30,pix-29,...pix+40+30]

 ;Prepare a big piece of spectrum for continuum fitting
   strr = reform(star(*,ordr))              ;template, one order
   flstrr = strr
   llo = max([0,place-180])
   hhi = min([llo+399,splen-1])
   llo = hhi-399
   dumsp = strr(llo:hhi)

  dumdum = dumsp
  lx = maxloc(dumdum(0:180),/first)
  dumdum(lx) = 0.
  lx = maxloc(dumdum(00:180),/first)
  rx = maxloc(dumdum(hhi-llo-180:hhi-llo),/last)+(hhi-llo-180)
  dumdum(rx) = 0.
  rx = maxloc(dumdum(hhi-llo-180:hhi-llo),/last)+(hhi-llo-180)
  cof = pl_fit([lx,rx],[dumsp(lx),dumsp(rx)],1)
  lcont = poly_fat(findgen(hhi-llo+1),cof)
  flstrr(llo:hhi) = dumsp/lcont

   str = reform(flstrr(lo:hi))     ;100 pixel chunk of flattened template
   wdum = lambda+dispersion*xxnow  ;rough wavelength scale

;Calculate the error in the chunk velocity analytically
;See EQN In Butler & Marcy 1996: Sigma(mean) = 1./SQRT(SUM(1/SIG^2))
  sp   = reform(normstar(place:place+len-1,ordr))
  eps  = sqrt(sp)
  didp = sp(1:len-1) - sp(0:len-2)       ;slope:   dI/d(pix)
  didv = didp*(lambda/(c*dispersion))    ;slope in real intensity per m/s
  dsst(n).weight = total((didv/eps)^2)   ;EQN 5 in error write up

;Making PSF with next line
   psfav,vd,ordr,place,osmp,psf,nip,del_ord=d_ord,del_pix=d_pix,orddist=orddst

;PB Kludge, Oct. 31, 1994  -- Bombing because no good PSF's are being found!
       if (nip lt 2) or (max(psf) le 0) then begin 
	  print,'DSST: problems finding good PSF!'
	  print,'DSST: increase averaging area, del_pix=100, del_ord=5'
	  psfav,vd,ordr,place,osmp,psf,nip,del_ord=d_ord+1,del_pix=d_pix+60
       endif

    if max(psf) le 0 then stop,'Max(PSF) < 0 in MDSST -- Flaw'
    if max(psf) eq min(psf) then stop,'Max(PSF) = Min(PSF) in MDSST -- Flaw'

    temp = str  ;template spectrum 
    wavtemp = lambda - dispersion*30. + findgen(100)*dispersion

    if m eq 0 then par(0)=3.d5   ;Bizarre Doppler shift: Flag 1st chunk 

;let dst be calculated for bad case, but give dsst.weight of -1, PB, 15Jun98
if n_elements(where(filt gt 0)) lt 20 then begin
       print,'Only '+strtrim(n_elements(where(filt gt 0)),2)+ $
	     ' good pixels in this chunk, DSST.WEIGHT = -1'
       dsst(n).weight=-1
       filt=intarr(n_elements(str))*0+1
endif

    deconv_m,par,plot=1      ;replaces deconv_j
;    deconv_j,str,psf,dum,filter=filt,bval=bval,parsm=4,osamp=osmp

    
    wv0  = lambda-dispersion*11.75  ;Allow for 11.75 extra pixels at edges.
;    wvlast = wv0 + 64.*dispersion   ;total length = 11.75 + 40 + 11.75 = 64
    i = where(wavfine ge wv0)
    wavfine = wavfine(i(0):i(0)+dstlen-1)     ;Force length = 256 or 296
    sunfine = sunfine(i(0):i(0)+dstlen-1)

;    wstr = findgen(n_elements(dum))*(dispersion/float(osmp)) + min(wdum)
;       px0=first_el(where(wstr ge wv0))              ;find pixel at 1st lambda
;       dum=dum(px0:px0+255)                          ;take 1st 256 osmp pixls
;       wstr=wstr(px0:px0+255)                        ;corresp. wavelengths
;       dsst(n).dst=dum                               ;deconv'd template!
;       dsst(n).w1=double(dispersion)/(double(osmp))

;Store in DSST Structure
dsst(n).dst = sunfine
dsst(n).w1 = wavfine(1)-wavfine(0)
dsst(n).w0 = wavfine(0)

;       if n_elements(wc) gt 1 then dsst(n).w0=wstr(0) else dsst(n).w0=wv0
    ;print,ordr,place,px0,lambda,wstr(0)

    IF keyword_set(plot) then begin
        xt='Wavelength ( '+ang()+' )'
        yt='Residual Intensity'
        ttl='Order '+strtrim(ordr,2)+', Pixel '+strtrim(place,2)
        plot,wavfine,sunfine,/xsty,/ynoz,xtitle=xt,ytitle=yt,title=ttl,co=135
        oplot,wavtemp,temp,co=180,ps=2
        oplot,wavfine,wavfine*0.+1.,co=191
          wait,1
    ENDIF
 IF n eq 0 then begin
   print,' '
   print,'                  *  CREATING DSST  *     '
   print,' '
   print,'|--------------------------------------------------------|'
   print,'| Order Pixel  Lambda   dlam/dx  Scat    SLOPE  # PSFs   |'
   print,'|               (A)     (A/pix)           SUM   avg''d    |'
   print,'|--------------------------------------------------------|'
 ENDIF
 fmt = '(I5,I6,F10.3,F9.5,F8.3,F8.2,I6)'
 print,format=fmt,ordr,place,lambda,dispersion,scat,dsst(n).weight,nip
			     
ENDFOR; n=0,(n_lines-1)
vd=oldvd

return
end  
