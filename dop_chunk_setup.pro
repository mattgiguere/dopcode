;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PRO DOP_CHUNK_SETUP
;
; CALLED FROM: DOP_MAIN
; CALLS TO: DOP_WAVOB
; CALLS TO: (OPTIONAL: DOP_FITWAV)
;
; PURPOSE:  
;   Set up an initial database (chunk_arr) structure
;   that will hold:
;     - order and pixel of chunk
;     - SOBS for the chunk
;     - SIOD and WIOD for the chunk
;     - ISS and WIS for the chunk
;     - SMOD: model spectrum
;     - FREE_PARs (wav soln, PSF, wav shift)  
;          WCOF: fltarr(2) = lambda and dispersion
;          Z: doppler shift
;          OFFSET: continuum shift
;          AMP[0]: misnomer - really the FWHM of the central Gaussian
;          AMP[1:16]: amplitudes of the flanking Gaussians (or nodes)
;     - FREE_IND boolean flags identifying free pars
;     - FIT: chisq fit between SOBS and SMOD
;     - IP: the PSF in oversampled pixels 
;     - WEIGHT: chunk weight from ISS - not from Doppler analysis!
;     - GDPIX: number of good pixels 
;
; PROCEDURE:
;   1) define chunk structure tag names for vdiodines
;   2) install wavelength soln from nearest iodine observation
;   3) adopt PSF from nearest iodine
;   4) initial doppler shift is the barycentric velocity of 
;      observation relative to the bary shift of the template 
;   5) (if dop_fitwav called) enforce continuity in the wavelengths 
;      across the order and update the initial wavelength coefficients 
;      obtained from the iodine observation 
;   6) (if dop_fitwav called) use this smoothed wavelength solution
;      for the initial wavelength solution
;   7) find the "pixob," where the first pixel of each 
;      ISS chunk corresponds in wavelength to the observation
;   8) return the loaded chunk structure to dop_driver
;
;  THREE POSSIBLE CASES (EACH HARDWIRED IN THIS PROGRAM): 
;   1) Program observations with dsst-style ISS
;   2) Iodine observations with FTS I2 
;   3) NSO Atlas for the ISS (daysky, moon, etc)
;
; INPUTS: 
;   OBSNM: name of iodspec to be analyzed
;   ISS_NM: name of deconvolved stellar template structure
;   ISS_BC: barycentric correction of ISS 
;   POLY_ORD: order of polynomial fit across order 
;   OSAMP: oversampling of pixels for PSF fitting
;   NSO: NSO spectrum as template 
;   DOPENV: structure with file paths, psf description
;   DCONV_SPEC: deconvolved spectrum with orders described 
;               by (wav, spec) array, like NSO
;   PLOT: optional plots
;   DEMO: more extensive optional plots
;   FITWAV: optional keyword to invoke fitting continuous 
;           wavelength across orders and to use that fit 
;           for the wavelength assignment in the chunk
; OUTPUTS:
;   CHUNK_ARR: returned chunk array, stocked with initial wavelength 
;          estimates and PSF values 
;   OBS_INFO: psfpix, psfsig, bccor, npix, issnm, iodnm
;
; OUTSTANDING ISSUES:
;   Need to make the ipcf.jd double precision
;   LATER: Expand the wavelength range to include orders 35 - 53?
;   Test passes: frz psf, frz disp and w0, check chunk gap sizes
;
; Written by: Debra Fischer, SFSU, Nov 2007
;             DF Jan 2013 implement median yrmo SLSF and disp for fiber-fed obs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO DOP_CHUNK_SETUP, dopenv, chunk_arr, obs_info, tag=tag, vdtag=vdtag, avg=avg,$
					 vdavgnm=vdavgnm, demo=demo, verbose=verbose, crosscorl=crosscorl, $
                     iodarr=iodarr, tmpl_dir=tmpl_dir, yrmo=yrmo,cdnear_name=cdnear_name

 ; DEFINE CONSTANTS AND VARIABLE NAMES
   scat_light=0.0  ;Fischer, for testing - not currently used
   c_light = dopenv.c_light
   osamp=dopenv.osamp
   if dopenv.iss_obnm eq 'nso' then nso = 1 else nso = 0
   iss_bc = dopenv.iss_bc
   xcorl_lambda = dopenv.xcorl_lambda
   ipcf_file = dopenv.ipcf_file
   baryfile = dopenv.bary_file       ; e.g., /mir1/bary/bcvel.dat
   restore, baryfile                 ; bcat structure
   restore, ipcf_file                ; ipcf structure
   zp=0                              ; zeroth pass doppler par's 
   fp=1                              ; first pass doppler par's
   lp=2                              ; last pass doppler par's
   npsf=n_elements(dopenv.psfpix)
   psfmod=dopenv.psfmod
   npix = dopenv.n_pix

if keyword_set(cdnear_name) then begin
	iodnm=cdnear_name
	p1=strpos(iodnm,'_')
	obnm=strmid(iodnm,p1+1,strlen(iodnm)-p1)
	x=where(bcat.obsnm eq obnm,nx) 
	if nx eq 0 then x=where(bcat.obsnm eq strmid(obnm,1,strlen(obnm)-1),nx)
	print, '====> Starname: ', dopenv.obj_nm
    print, '====> Observation: ',dopenv.obsnm
endif

if ~keyword_set(cdnear_name) then begin
 ; FIND THE CLOSEST (IN TIME) IODINE WAVELENGTH AND PSF SOLN
 ; AND RETURN THAT IODNM TO DOP_DRIVER   
   x=where(bcat.obsnm eq 'r'+strt(strmid(dopenv.obsnm, 1, strlen(dopenv.obsnm))),nx) 
   if nx le 0 then begin
;   x=where(strt(bcat.obsnm) eq strt(strmid(dopenv.obsnm, 0, strlen(dopenv.obsnm))),nx) 
   if strmid(dopenv.obsnm,0,3) eq 'chi' then x=where(strt(bcat.obsnm) eq dopenv.obsnm,nx) 
   if strmid(dopenv.obsnm,0,4) eq 'achi' then x=where(strt(bcat.obsnm) eq strmid(dopenv.obsnm,1,strlen(dopenv.obsnm)-1),nx)
   if nx eq 0 then stop 
   endif
                    
   print, '====> Starname: ', bcat[x].objnm
   print, '====> Observation: ',dopenv.obsnm

   diff=double((ipcf.jd-2440000.)-bcat[x].jd)  ; diff is an array of delta(JD) 
   adiff=abs(diff) 
   n_ind = where(adiff eq min(adiff))  &  n_ind = n_ind[0]
   ;vdiodnm='vd'+tag+ipcf[n_ind].iodnm+'_'+ipcf[n_ind].obnm 
   iodnm='cd'+vdtag+'b'+ipcf[n_ind].iodnm+'_'+ipcf[n_ind].obnm 
endif ;find closest in time 

   fvd=findfile(dopenv.files_dir+iodnm,count=count_iodnm)
       if count_iodnm eq 0 and ~keyword_set(cdnear_name) then begin
       		print,'running dop_create_ipcf...'
       		dop_create_ipcf,tag=tag,observatory=dopenv.observatory
       		restore,ipcf_file
   			diff=double((ipcf.jd-2440000.)-bcat[x].jd)  ; diff is an array of delta(JD) 
   			adiff=abs(diff) 
   			n_ind = where(adiff eq min(adiff))  &  n_ind = n_ind[0]
   			;vdiodnm='vd'+tag+ipcf[n_ind].iodnm+'_'+ipcf[n_ind].obnm 
   			iodnm='cd'+vdtag+'b'+ipcf[n_ind].iodnm+'_'+ipcf[n_ind].obnm 
       		fvd=findfile(dopenv.files_dir+iodnm,count=count_iodnm)
       		if count_iodnm eq 0 then stop, 'vdiod file is missing'
       	endif

   print, '====> Closest match cdbiod: '+iodnm
   print, '  '

   restore,fvd      ; iodine soln for initial wavelength and PSF 
   vd_near = chunk_arr
     
   xord=vd_near.ordt
   gg=uniq(xord)
   ordlist=xord[gg]
   
   ; LIGHT SMOOTHING / CONTINUITY OF WAVELENGTH AND DISPERSION FOR ZP 
   for k=0,n_elements(ordlist)-1 do begin
        xord=where(vd_near.ordob eq ordlist[k],nx)
                sub_chunk=vd_near[xord]
                ; smooth the lp dispersion and wavelength for the initial guess
                dispsm=dop_smooth(ordlist[k], sub_chunk, pass=lp, /disp) 
                wavsm=dop_smooth(ordlist[k], sub_chunk, pass=lp, /wav) 
                for ii=0, nx-1 do (*vd_near[xord[ii]].free_par)[zp].wcof[1] = dispsm[ii]
                for ii=0, nx-1 do (*vd_near[xord[ii]].free_par)[zp].wcof[0] = wavsm[ii]
                for ii=0, nx-1 do (*vd_near[xord[ii]].free_par)[fp].wcof[1] = dispsm[ii]
                for ii=0, nx-1 do (*vd_near[xord[ii]].free_par)[fp].wcof[0] = wavsm[ii]
   endfor
                

 ; REPLACE PSF PARS WITH MEDIANS FROM THE YRMO SET OF BSTARS
   if keyword_set(avg) then begin
   ;fischer jan7, 2013
   	  if ~keyword_set(vdavgnm) then vd_avgfile='vd'+tag+'AVG_'+yrmo else vd_avgfile=vdavgnm
   	  ;vd_avgfile='vdargAVGs_1212'
   	  fvdav=file_search(dopenv.files_dir+vd_avgfile,count=ccount)
   	  if ccount eq 0 then stop,'no vd_avgfile found'
   	  restore,dopenv.files_dir+vd_avgfile  ;chunk_avg
   	  if psfmod eq 'gaussian' then begin
   	  	for ij=0,dopenv.n_chunks-1 do begin
   	  		vd_near[ij].ip[*,1]=vdav[ij].psf  ;chunk_avg[ij].ip[*,1] first pass
   	  		vd_near[ij].ip[*,2]=vdav[ij].psf  ;chunk_avg[ij].ip[*,2] second pass
   	  		for jk=0,npsf-1 do (*vd_near[ij].free_par)[1].amp[jk]=vdav[ij].par[jk] 
   	  		for jk=0,npsf-1 do (*vd_near[ij].free_par)[2].amp[jk]=vdav[ij].par[jk] 
	  	endfor
	  endif
;	  if psfmod eq 'bspline' then begin 
;	  endif 	
   endif
	
   chunk_arr=vd_near
   	
   if keyword_set(verbose) then begin
      print,(*vd_near[0].free_par).wcof[*]
      print,(*vd_near[0].free_par).amp[*]
   endif

   count_chunk=0
   nch=n_elements(vd_near.ordt)
   xcorl_disp=(*vd_near[nch-1].free_par)[lp].wcof[1]  ;close enough?

 ; OBS_INF0 CONTAINS INFORMATION THAT IS OBSERVATION-SPECIFIC
 ; RATHER THAN CHUNK-SPECIFIC
   	if psfmod eq 'gaussian' then begin
   		obs_info={psfpix:dopenv.psfpix, $ ;PSF parameters: positions of Gaussians
             psfsig:dopenv.psfsig, $ ;PSF parameters: FWHM of Gaussians
             bccor:bcat[x].bc,     $ ;barycentric correction used in calculating vel
             npix:npix,    $ ;e.g.,40 or 80-pixel chunks
             issnm:dopenv.iss_nm,  $
             decker:dopenv.decker, $ ;e.g. B1
             iodnm:iodnm}            ;iodine obs used for initial par guesses
	endif
    if psfmod eq 'bspline' then begin
   		obs_info={psfbsplnplaces:dopenv.psfbsplnplaces, $ ;PSF parameters: positions 
   			 psfbsord:dopenv.psfbsord,  $ ;bspline order 
             bccor:bcat[x].bc,     $ ;barycentric correction used in calculating vel
             npix:npix,    $ ;e.g.,40 or 80-pixel chunks
             issnm:dopenv.iss_nm,  $
             decker:dopenv.decker, $ ;e.g. B1
             iodnm:iodnm}           ;iodine obs used for initial par guesses       
	endif
             

 ; BEGIN CASE 1 **************************************************************
 ; PROGRAM OBSERVATIONS WITH DSST-STYLE ISS
   if dopenv.iss_nm ne 'iod' and dopenv.iss_nm ne 'nso' then begin
      restore, dopenv.files_dir+dopenv.iss_nm      
         test=where(tag_names(dsst) eq 'W0',ntest)  ;dsst instead of iss? 
         if ntest gt 0 then iss=dsst
      nchunks=n_elements(iss.ordt)
      chunk_arr.weight = iss.weight
      if (iss[1].pixt - iss[0].pixt) ne npix then stop,'Mismatched DSST'
      is_spec=mrdfits(tmpl_dir+dopenv.iss_obnm+'.fits',dopenv.file_ext)
      if n_elements(size(is_spec)) gt 5 then is_spec=reform(is_spec[1,*,*])  ;wavelength array was attached
      ob_spec=mrdfits(dopenv.obs_dir+dopenv.obsnm+'.fits',dopenv.file_ext)
      if n_elements(size(ob_spec)) gt 5 then ob_spec=reform(ob_spec[1,*,*])  ;wavelength array was attached
    ; start with BC velocity correction relative to the ISS observation
      init_z = double((iss_bc - obs_info.bccor)/c_light) 

    ; QUICK XCORL OUTSIDE THE I2 REGION TO REFINE THE INITIAL WAVELENGTH SHIFT
    ; BETWEEN THE TEMPLATE AND OBSERVATION 
      if keyword_set(crosscorl) then begin 
      ; THIS LOOP IS NEEDED FOR SB'S - OTHERWISE USE INIT_Z
          snip_iss=is_spec[dopenv.xcorl_pix1:dopenv.xcorl_pix1+dopenv.xcorl_npix,dopenv.xcorl_ord]
          snip_obs=ob_spec[dopenv.xcorl_pix1:dopenv.xcorl_pix1+dopenv.xcorl_npix,dopenv.xcorl_ord]
          contf,snip_iss,c_iss,sbin=20,nord=2
          contf,snip_obs,c_obs,sbin=20,nord=2
          shift_range=dopenv.xcorl_shiftrange ; pixel range for xcorl 
          xcorlb,snip_obs/c_obs, snip_iss/c_iss, shift_range, pix_shift
          init_z = double(pix_shift*xcorl_disp/xcorl_lambda)
          print,init_z
          if keyword_set(demo) then begin
             plot,snip_obs/c_obs,col=1,title='Xcorl: template and observation'
             ;oplot,snip_iss/c_iss, col=155  ; unshifted
             shfour,snip_iss/c_iss,pix_shift,new_iss
             oplot,new_iss,col=222
             print,'shift ISS by : ',pix_shift
          endif
       endif

      for i=0,nchunks-1 do begin  
       ; first guess model should inherit ISS wavelengths (will be shifted I2 wavelengths), 
       ; dispersion, and the VDIOD psf
       ; the dsst must have been generated with the same dopenv.dpad used here.
       ;;;;;;fischer jan 13, 2013 change following 2 lines (use bstar wave soln)
         ;(*chunk_arr[i].free_par)[zp].wcof[0] = iss[i].w0 + dopenv.dpad*osamp*iss[i].w1
         ;(*chunk_arr[i].free_par)[zp].wcof[1] = iss[i].w1*osamp    
		 ;(*chunk_arr[i].free_par)[zp].wcof[0] = (*vd_near[i].free_par)[lp].wcof[0]
 		 ;(*chunk_arr[i].free_par)[zp].wcof[1] = (*vd_near[i].free_par)[lp].wcof[1]
 		
         if dopenv.n_wcof eq 3 then (*chunk_arr[i].free_par)[zp].wcof[2] = iss[i].w2*osamp 
       ; only the first pass has amplitudes - the last pass uses smoothing 
         ;for j=0, npsf-1 do (*chunk_arr[i].free_par)[zp].amp[j] = (*vd_near[i].free_par)[fp].amp[j]  
         (*chunk_arr[i].free_par)[zp].z = init_z
         (*chunk_arr[i].free_par)[fp].z = init_z
      endfor 

    ; MATCH WAVELENGTH CHUNKS BETWEEN THE TEMPLATE AND PROGRAM OBSERVATION
    ; ISS_DISP * (PIXT - PIXOB) = DELTA(LAMBDA) = Z * ISS_W0
    ; PIXOB = PIXT - Z * (ISS_W0 / ISS_DISP)
      ordr = chunk_arr.ordob              
      gg = uniq(ordr) 
      ordr = ordr(gg)            ; array of orders in the chunk
      nord = n_elements(ordr)    ; number of orders
	  
    ; ONE ORDER AT A TIME...
      for i=0,nord - 1 do begin 
         obs_ord = reform(ob_spec[*,ordr[i]]) 
         ch_ord = where(chunk_arr.ordob eq ordr[i],nch_ord)   
         ncol=n_elements(ob_spec[*,ordr[i]])
         subch = chunk_arr[ch_ord]              ; with iss wavelengths 
         
;       ; CHECK FOR ANOMALOUS PSF
;         ck_chunkpsf,subch_in=subch,subch_out=subch_out
;         subch=subch_out  ; replace or not               
		 tmp=subch 
		 
;;;;;THIS COULD BE A PROBLEM? - forces wavelength continuity   
          wav=dop_wavob(tmp, ncol,poly_ord=3)  ; wavelength for every pixel
		 
;		 wav=fltarr(ncol)
;		 ;pixels 0 - 79
;		 wfirst=double((*subch[0].free_par)[zp].wcof[0] - npix*(*subch[0].free_par)[zp].wcof[1])
;		 wav[0:npix-1] = wfirst + double(indgen(npix)*(*subch[0].free_par)[zp].wcof[1])
;		 for kk = 0, nch_ord-1 do begin  ; pixels 79:3020
;		 	ll=indgen(npix)+((kk+1)*npix)
;		 	; using the dop_smoothed input for each independent chunk
;		 	wav[ll]=double((*subch[kk].free_par)[zp].wcof[0] + indgen(npix)*(*subch[kk].free_par)[zp].wcof[1])
;		 endfor ;kk
;		 ;pixels 3140 - 3200
;		 wlast=max(wav[ll]) 
;		 ind1=max(ll)+1
;		 wav[ind1:ind1+79]=wlast + (1+indgen(npix))*(*subch[nch_ord-1].free_par)[zp].wcof[1]
;stop
       ; FIND WAVELENGTHS FOR EACH PIXEL IN THE ORDER
         for kk=0,nch_ord-1 do begin
          ; GET PIXOB (BARY-SHIFTED or DOPPLER-SHIFTED)
            shft_wav = double((*subch[kk].free_par)[zp].wcof[0]) + init_z*double((*subch[kk].free_par)[zp].wcof[0])
            if dopenv.n_wcof eq 3 then shft_wav = double((*subch[kk].free_par)[zp].wcof[0]) + $
            	init_z*double((*subch[kk].free_par)[zp].wcof[0])

            shft_wav = double((*subch[kk].free_par)[fp].wcof[0]) + init_z*double((*subch[kk].free_par)[fp].wcof[0])
            if dopenv.n_wcof eq 3 then shft_wav = double((*subch[kk].free_par)[fp].wcof[0]) + $
            	init_z*double((*subch[kk].free_par)[fp].wcof[0])

            diff = abs(shft_wav - wav)
            pix_ind = where(diff eq min(diff))   ;pixel with closest wavelength 
            subch[kk].pixob = pix_ind          	
          	  
          ; RESET W0 FOR EACH CHUNK (PASS = 0), BASED ON PIXOB
            start_pix = subch[kk].pixob
            pix1=start_pix  &  pix2=start_pix+npix-1

          ; FILL (*CHUNK.SOBS) WITH OBSERVED SPECTRUM
            (*subch[kk].sobs) = obs_ord[pix1:pix2]
            (*subch[kk].free_par)[zp].wcof[0] = double( wav[pix1] )
            (*subch[kk].free_par)[zp].wcof[1] = double( (wav[pix2]-wav[pix1])/npix )
            (*subch[kk].free_par)[fp].wcof[0] = double( wav[pix1] )
            (*subch[kk].free_par)[fp].wcof[1] = double( (wav[pix2]-wav[pix1])/npix )
            subch[kk].cts = long(median((*subch[kk].sobs)))

          ; FILL (*CHUNK.SIOD)
            ipad=1.2  ; extra wavelength
            read_fts,wav[pix1]-ipad,wav[pix2]+ipad, dopenv, wiod,siod
;            siod = siod - (scat_light*siod)
            contf,siod,ciod,sbin=20,nord=1
            (*subch[kk].siod) = siod/ciod
            (*subch[kk].wiod) = wiod

          ; FILL (*CHUNK.SIS)
            iss_chunk=count_chunk
            (*subch[kk].sis) = iss[iss_chunk].dst
            darr = dindgen(n_elements((*subch[kk].sis)))
            (*subch[kk].wis) = iss[iss_chunk].w0 + iss[iss_chunk].w1*darr
            count_chunk=count_chunk+1

         endfor  ; kk: all chunks in this order
         chunk_arr[ch_ord] = subch  ; update chunk_arr
      endfor    ; i: each order

   endif  ; program star observation

 ; END CASE 1 **************************************************************
 ; PROGRAM OBSERVATIONS WITH DSST-STYLE ISS

 ; BEGIN CASE 2 **************************************************************
 ; IODINE OBSERVATIONS WITH FTS-I2 FOR THE ISS
   if dopenv.iss_nm eq 'iod' then begin
      nchunks=n_elements(vd_near.ordt)
      chunk_arr = dop_chunk_replicate(nchunks, dopenv)
      chunk_arr.ordt = vd_near.ordt
      chunk_arr.ordob = vd_near.ordt  ; no shift for iodine
      chunk_arr.pixt = vd_near.pixt
      chunk_arr.pixob = vd_near.pixob
      chunk_arr.weight = 1.0
      init_z=0.0d

    ; INITIALIZE THE ZEROTH [0] ITERATION 
      for i=0,nchunks-1 do begin  
         (*chunk_arr[i].free_par)[zp].wcof[0] = (*vd_near[i].free_par)[lp].wcof[0]
         (*chunk_arr[i].free_par)[zp].wcof[1] = (*vd_near[i].free_par)[lp].wcof[1]
         (*chunk_arr[i].free_par)[zp].amp[0] = dopenv.psfsig[0]
         for j=1, 10 do (*chunk_arr[i].free_par)[zp].amp[j] = (*vd_near[i].free_par)[fp].amp[j]
         (*chunk_arr[i].free_par)[zp].z = init_z

         (*chunk_arr[i].free_par)[fp].wcof[0] = (*vd_near[i].free_par)[lp].wcof[0]
         (*chunk_arr[i].free_par)[fp].wcof[1] = (*vd_near[i].free_par)[lp].wcof[1]
         (*chunk_arr[i].free_par)[fp].amp[0] = dopenv.psfsig[0]
         for j=1, 10 do (*chunk_arr[i].free_par)[fp].amp[j] = (*vd_near[i].free_par)[fp].amp[j]
         (*chunk_arr[i].free_par)[fp].z = init_z
      endfor

      ordr = chunk_arr.ordob              
      gg = uniq(ordr) 
      ordr = ordr(gg)            ; array of orders in the chunk
      nord = n_elements(ordr)    ; number of orders
      ob_spec=mrdfits(dopenv.obs_dir+dopenv.obsnm+'.fits',dopenv.file_ext)
	  if n_elements(size(ob_spec)) gt 5 then ob_spec=reform(ob_spec[1,*,*])  ;wavelength array was attached

      ; ONE ORDER AT A TIME...
      for i=0,nord - 1 do begin 
         ch_ord = where(chunk_arr.ordob eq ordr[i],nch_ord)   
         subch = chunk_arr[ch_ord]     ; one order at a time
         ncol=n_elements(ob_spec[*,ordr[i]])
         
        ; CHECK FOR ANOMALOUS PSF
       ;  ck_chunkpsf,subch_in=subch,subch_out=subch_out
       ;  subch=subch_out  ; replace or not               
		 tmp=subch      
		 
;;;;THIS COULD BE A PROBLEM
         wav=dop_wavob(subch, ncol, poly_ord=6.)  ; one wavelength for every pixel

;		 wav=fltarr(ncol)
;		 ;pixels 0 - 79
;		 wfirst=double((*subch[0].free_par)[zp].wcof[0] - npix*(*subch[0].free_par)[zp].wcof[1])
;		 wav[0:npix-1] = wfirst + double(indgen(npix)*(*subch[0].free_par)[zp].wcof[1])
;		 for kk = 0, nch_ord-1 do begin  ; pixels 79:3020
;		 	ll=indgen(npix)+((kk+1)*npix)
;		 	; using the dop_smoothed input for each independent chunk
;		 	wav[ll]=double((*subch[kk].free_par)[zp].wcof[0] + indgen(npix)*(*subch[kk].free_par)[zp].wcof[1])
;		 endfor ;kk
;		 ;pixels 3140 - 3200
;		 wlast=max(wav[ll]) 
;		 ind1=max(ll)+1
;		 wav[ind1:ind1+79]=wlast + (1+indgen(npix))*(*subch[nch_ord-1].free_par)[zp].wcof[1]

         obs_ord = reform(ob_spec[*,ordr[i]]) 

         for kk=0,nch_ord-1 do begin
            start_pix = subch[kk].pixob
            pix1=start_pix  &  pix2=start_pix+npix-1

          ; FILL (*CHUNK.SOBS) WITH OBSERVED SPECTRUM
            (*subch[kk].sobs) = obs_ord[pix1:pix2]
            subch[kk].cts = long(median((*subch[kk].sobs)))

          ; FILL (*CHUNK.SIOD)
            ipad=1.2  ; extra wavelength
            read_fts,wav[pix1]-ipad,wav[pix2]+ipad, dopenv, wiod,siod
;            siod = siod - (scat_light*siod)
            contf,siod,ciod,sbin=20,nord=1
            (*subch[kk].siod) = siod/ciod
            (*subch[kk].wiod) = wiod
            narr=n_elements(siod)

          ; FILL (*CHUNK.SIS) = 1.0 ARRAY
            (*subch[kk].sis) = fltarr(narr)*0.0d + 1.0d
            (*subch[kk].wis) = wiod

         endfor  ; kk: all chunks in this order
         chunk_arr[ch_ord] = subch 
      endfor
   endif

 ; END CASE 2 **************************************************************
 ; IODINE OBSERVATIONS WITH FTS-I2 FOR THE ISS

 ; BEGIN CASE 3 **************************************************************
 ; OBSERVATIONS WITH NSO ATLAS FOR THE ISS (MOON, DAYSKY, ETC)
   ;;;;THIS BLOCK OF CODE NEEDS WORK
   if dopenv.iss_nm eq 'nso' then begin
      nchunks=n_elements(vd_near.ordt)
      chunk_arr = dop_chunk_replicate(nchunks, dopenv)
      chunk_arr.ordt = vd_near.ordt
      chunk_arr.ordob = vd_near.ordt  ; no shift for iodine
      chunk_arr.pixt = vd_near.pixt
      chunk_arr.pixob = vd_near.pixob
      chunk_arr.weight = 1.0
;      init_z = bcat[x].bc / c_light

    ; QUICK XCORL TO FIND THE INITIAL WAVELENGTH SHIFT
    ; BETWEEN THE TEMPLATE AND OBSERVATION 
          ob_spec=mrdfits(dopenv.obs_dir+dopenv.obsnm+'.fits',dopenv.file_ext)
          if n_elements(size(ob_spec)) gt 5 then ob_spec=reform(ob_spec[1,*,*])  ;wavelength array was attached
          ; Begin CTIO
          ;;; I know this is bad code.... 
          snip_obs=ob_spec[1000:1300,3]  ; 4967 - 4984 A
          restore,dopenv.files_dir+'ctio_rqa06.7400.dat'
          snip_w=w[1000:1300,3]  
          rdnso,w_nso,s_nso, 4966, 4985
          par=2.0
          ip=psf(par)
          con_nso=convol(s_nso, ip, /edge_truncate, /normalize)
          rebin, w_nso, con_nso, snip_w, snip_iss
          ; end CTIO kludge
          ; Lick
;          snip_obs=ob_spec[1000:1300,55]
;          rdsi,is_spec,dopenv.iss_obnm
;          snip_iss=is_spec[1000:1300,55]
;          contf,snip_iss,c_iss,sbin=20,nord=2
          
;          contf,snip_obs,c_obs,sbin=20,nord=2
;          shift_range=25        ; pixel range for xcorl 
;          xcorlb,snip_obs/c_obs, snip_iss/c_iss, shift_range,
;          pix_shift
;          xcorlb,snip_obs/c_obs, snip_iss, shift_range, pix_shift
;          plot,snip_w,snip_obs/c_obs,/xsty, col=1
;          oplot,snip_w,snip_iss, col=222
;          pix_shift=0.
;          print,pix_shift
;          init_z = double(pix_shift*xcorl_disp/xcorl_lambda)
          init_z = 0.

    ; INITIALIZE THE ZEROTH [0] ITERATION 
      for i=0,nchunks-1 do begin  
         (*chunk_arr[i].free_par)[zp].wcof[0] = (*vd_near[i].free_par)[lp].wcof[0]
         (*chunk_arr[i].free_par)[zp].wcof[1] = (*vd_near[i].free_par)[lp].wcof[1]
         (*chunk_arr[i].free_par)[zp].amp[0] = dopenv.psfsig[0]
         for j=1, 10 do (*chunk_arr[i].free_par)[zp].amp[j] = (*vd_near[i].free_par)[fp].amp[j]
         (*chunk_arr[i].free_par)[zp].z = init_z
      endfor

      ordr = chunk_arr.ordob              
      gg = uniq(ordr) 
      ordr = ordr(gg)            ; array of orders in the chunk
      nord = n_elements(ordr)    ; number of orders

      ; ONE ORDER AT A TIME...
      for i=0,nord - 1 do begin 
         ch_ord = where(chunk_arr.ordob eq ordr[i],nch_ord)   
         subch = chunk_arr[ch_ord]     ; one order at a time
         ob_spec=mrdfits(dopenv.obs_dir+dopenv.obsnm+'.fits',dopenv.file_ext)
         if n_elements(size(ob_spec)) gt 5 then ob_spec=reform(ob_spec[1,*,*])  ;wavelength array was attached
;         rdsi,ob_spec,dopenv.obsnm
         ncol=n_elements(ob_spec[*,ordr[i]])
         wav=dop_wavob(subch, ncol, poly_ord=6)  
         obs_ord = reform(ob_spec[*,ordr[i]]) 
 
         for kk=0,nch_ord-1 do begin
            start_pix = subch[kk].pixob
            pix1=start_pix  &  pix2=start_pix+npix-1

          ; FILL (*CHUNK.SOBS) WITH OBSERVED SPECTRUM
            (*subch[kk].sobs) = obs_ord[pix1:pix2]
            subch[kk].cts = long(median((*subch[kk].sobs)))

          ; FILL (*CHUNK.SIOD)
          ipad=1.2
          read_fts,wav[pix1]-ipad,wav[pix2]+ipad, dopenv, wiod,siod
;            siod = siod - (scat_light*siod)
            contf,siod,ciod,sbin=20,nord=1
            (*subch[kk].siod) = siod/ciod
            (*subch[kk].wiod) = wiod
            narr=n_elements(siod)

          ; FILL (*CHUNK.SIS) = NSO
            vactoair,wmin  & vactoair,wmax
            rdnso, wis, sis, wmin, wmax
            airtovac, wis 
            (*subch[kk].sis) = sis
            (*subch[kk].wis) = wis
         endfor  ; kk: all chunks in this order
         chunk_arr[ch_ord] = subch 
      endfor
   endif

 ; END CASE 3 **************************************************************
 ; OBSERVATIONS WITH NSO ATLAS FOR THE ISS (MOON, DAYSKY, ETC)

   chunk_arr.gdpix = npix


; don't do this
;   if keyword_set(fitwav) then begin
;      ; plot = 1 (RMS scatter for each order)
;      ; plot = 2 (left and central chunk points along each order
;      ; demo = 1 (as an example, shows fitting for first order)
;      print,'   '
;      print,'====>Be patient - fitting for all orders'
;      print,'   '
;      new_chunk = dop_fitwav(chunk_arr, poly_ord=poly_ord, pass=0, plot=0, demo=0)
;      chunk_arr = new_chunk
;   endif
      
 ; THE STOCKED VD CAN BE RETURNED TO DOP_DRIVER NOW  


end


