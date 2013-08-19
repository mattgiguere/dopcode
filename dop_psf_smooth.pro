;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_PSF_SMOOTH
;
; initial coding: Debra Fischer, Feb 2008
;
; CALLED FROM: DOP_FIT2 
;       restore,'/mir1/files/vdaHIP8102_rg92.234' 
;       chunk=chunk_arr
;       dopenv=dop_init('rg92.234', 'dsstHIP8102q_rg92.dat',4494.078)
;       ip_av=dop_psf_smooth(chunk, 42, 850, dopenv=dopenv)
;
; PURPOSE: 
;   This function returns an averaged, smooth PSF based on the first pass
;   Doppler analysis using: 
;        psfpar[4] is the width of the central Gaussian for the PSF
;        psfpar[5] - psfpix[14] are amplitudes of little Gaussians for the PSF
;
;   Other free pars:
;        psfpar[0] is w0, the wavelength of the first pixel in the chunk
;        psfpar[1] is the linear dispersion coefficient
;        psfpar[2] is the Doppler shift, z = v/c
;        psfpar[3] is a scaling term for the continuum
;
; PROCEDURE
;   1) restore the first-pass PSF amplitudes stored in the output vd
;   2) calculate a median LSF for chunks with good fits.  
;
; INPUTS:
;   CHUNK: velocity data structure 
;   ORDER: order for smoothing
;   PIXEL: pixel for smoothing
;   OSAMP: the pixel oversample
;   DOPENV: a structure with
;           del_ord: the range of orders for averaging
;           del_pix: the range of pixels for averaging
;           psfpix: the position of Gaussians
;           psfsig: the widths of Gaussians
;           osamp: oversampling used with the DST
;
; OUTPUTS:
;   IP: instrumental profile 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_PSF_SMOOTH, chunk, order, pixel, pass=pass, $ 
                    dopenv=dopenv, chunk_arr=chunk_arr, plot=plot

  del_pix = dopenv.del_pix
  del_ord = dopenv.del_ord
  sigpsf = dopenv.psfsig
  pixpsf = dopenv.psfpix
  npix = dopenv.n_pix
  ord_dist = dopenv.ord_dist     ; approx num pixels between orders
  osamp = dopenv.osamp
  fp = 1   ; first Doppler pass par's - always for dop_psf_smooth

  xsm = where(chunk_arr.ordob le (order + del_ord) and $
              chunk_arr.ordob ge (order - del_ord) and $
              chunk_arr.pixob le (pixel + del_pix) and $
              chunk_arr.pixob ge (pixel - del_pix),nxsm)

  if nxsm eq 0 then stop, 'no chunks for PSF_smooth'

  chunk_sm = chunk_arr[xsm]
  chunk_ip=chunk.ip[*,1]

;  if ~keyword_set(pass) then pass=0
;  if pass eq 3 then begin  ;no first pass fit, so kludge this
;  	chunk_sm.fit[1] = chunk.fit[0]
;  endif

  ; CHUNK WEIGHTS FOR PSF AVERAGING
  ord_sep = abs(chunk_sm.ordob - order)*ord_dist
  ord_wt = 1. / sqrt(1.0*ord_sep/npix + 1.)  ; weight per chunk
  pix_sep = abs(chunk_sm.pixob - pixel)
  pix_wt = 1. / sqrt(pix_sep/npix + 1.)

  fit_wt=fltarr(nxsm)
  for i=0,nxsm-1 do fit_wt[i]=1./chunk_sm[i].fit[1]

  wt = ord_wt * pix_wt * fit_wt

  ; ONE IP FOR EACH CHUNK IN AVERAGING REGION
  num=30                        ; num*osamp = range for psf
  xarr=findgen(num*osamp + 1)/ osamp - num/2
  nkernal=n_elements(xarr)
  ip_array=fltarr(nkernal,nxsm)    ; e.g., ip(121,9)
  ip_av=fltarr(nkernal)

  ; WEIGHT EACH OF THE CHUNK PSF'S 
  ; NOTE!!  THE "FREE" IP (TO BE AVERAGED) IS HARDWIRED 
  ;         AS [*,0] FOR THE FIRST PASS
for i=0,nxsm-1 do begin
  ip_one=chunk_sm[i].ip[*,1]  ;one of the nearest neighboring chunk psfs

 ; if dopenv.psfmod eq 'gaussian' then begin
  ;     ; SHIFT TO CENTER - CHECK THIS - MAY WANT TO ELIMINATE!
;  fwhm=0.5*max(ip_one)
;  xx=where(ip_one ge fwhm and abs(xarr) lt 5., nxx)   ; peak points
 ; xx=where(ip_one eq max(ip_one), nxx)   ; peak points
;  if nxx ge 3 then begin   ; shift to center
;	dd = where(ip_one[xx] eq max(ip_one[xx]))  &  dd=dd[0]
 ;  if nxx le 2 then cntr=xarr[xx[0]]  
 ;  if nxx eq 3 then cntr=xarr[xx[1]]
;	if abs(cntr) ge 0.1 and abs(cntr) lt 1.8 then begin
;	   if dopenv.psfmod eq 'gaussian' then begin
;		 ip_one = dop_psf((*chunk_sm[i].free_par)[fp].amp, dopenv=dopenv, xarr=xarr, cntr=cntr)
;	   endif;gaussian
	;   if dopenv.psfmod eq 'bspline' then begin
  ;            print, (*chunk_sm[i].free_par)[fp].amp[1:*]
 ;             ip_one = dop_psf_bspline((*chunk_sm[i].free_par)[fp].amp[1:*], dopenv=dopenv, xarr=xarr, cntr=dopenv.cntr)
;	   endif;bspline
;	endif; 0.1 < cntr < 1.2
;	if abs(cntr) ge 1.8 then begin  ;dump far-out IP's
;	   wt[i]=0.
;	   ip_one= chunk_ip
;	endif
;  endif  ; shift to center
;  endif; PSF centering

  ip_array[*,i]=ip_one*wt[i]
endfor

  for j = 0, nkernal-1 do ip_av[j] = median(ip_array[j,*])

  ip_av = osamp * ip_av / total(ip_av)
 
  if keyword_set(plot) then begin
     plot,xarr,ip_av, col=1,/xsty, xra=[-7,7]
  endif

return, ip_av

end


