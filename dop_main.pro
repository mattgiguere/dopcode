;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PRO DOP_MAIN
;
; PROCEDURES CALLED: DOP_INIT, DOP_CHUNK_SETUP, DOP_MARQ
;
; PURPOSE: 
;   Engine for Doppler analysis
;   Based on algorithms described in: 
;      Marcy & Butler, 1992 PASP 104 270
;      Butler et al. 1996 PASP 108, 500 
;      Valenti et al. 1995 PASP 107 966
;
;
; PROCEDURE:
;   1) set up the pathnames and input variables:
;      DOP_INIT takes program-specific info (e.g., CTIO_INIT.PRO)
;   2) set up initial velocity data structure (DOP_CHUNK_SETUP)
;   3) starting with the template or ISS and FTS atlas, synthesize a 
;      nightly observation with Levenberg-Marquardt search for 
;      best free parameters
;      - restore the observation
;      - restore the ISS 
;      - step through each order
;      - slice and dice each order into wavelength chunks
;      - set pixel weights and filter out bad pixels 
;      - extract a matching segment from the FTS and normalize
;      - for each chunk, shift, multiply the ISS by the FTS, 
;        and convolve with PSF until best fit model
;
; CALLING SEQUENCE: 
;   dop_main, '128620',dopenv, obsnm='rg92.234', pass=pass, $
;             tag=tag, lm=lm
;
; INPUTS: 
;   STARNAME 
;   DOPENV structure (from ctio_init and dop_init):
;       OBSNM:  (obs name) stellar spectrum to be analyzed
;       ISS_NM: name of deconvolved template
;           structure has keywords: 
;             ordt, pixt, w0, disp, weight, iss[256]
;       ISS_BC: barycentric velocity of the template observation
;           (used to estimate the initial Doppler shift)
;   PASS: 1 or 2 
;   TAG: for filename of output vd and cd structures
;   LM: the type of Levenberg-Marquardt fitter
;       "MPF" = mpfit programs from Craig Markwardt
;       "STARSOLVE" = old L-M fitter written by J. Valenti 
;
; OPTIONAL INPUT: 
;   NSO: NSO atlas 
;   DCONV_SPEC: deconvolved spectrum in the format (wav, spec) 
;               for each order - my hope for future ISS's
;
; OUTPUTS
;   CHUNK_ARR: a structure with model parameters for each
;     wavelength segment (a.k.a. chunk) of the observed spectrum
;     stored in the complete data base or "cdxa" and "cdxb" files
;   VD: velocity data structure, a "light-weight" version of 
;     the "cdxb" file
;
; Written by Debra Fischer, SFSU, Nov 2007
; John Johnson, IfH, May 2008 (implementation of IDL Pointers)
;
; OUTSTANDING: 
;   IS VEL UPDATED FOR SEQUENTIAL GUESSES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO DOP_MAIN, starname, dopenv, obsnm=obsnm, pass=pass, tag=tag, nso=nso, $
              vdtag=vdtag, pix_filt = pix_filt, smdisp=smdisp, smwav=smwav, $
              iod_soln=iod_soln, plot=plot, tmpl_dir=tmpl_dir, avg=avg, vdavgnm=vdavgnm,  $
              yrmo=yrmo, demo=demo, cdnear_name=cdnear_name, verbose=verbose, $
              lm=lm, crosscorl=crosscorl

 iss_nm=dopenv.iss_nm
 iss_bc=dopenv.iss_bc
 iss_obnm=dopenv.iss_obnm
 if ~keyword_set(smdisp) then smdisp=0  
 if ~keyword_set(smwav) then smwav=0
 if ~keyword_set(vdtag) then vdtag=tag
 
 zp = 0                      ; zeroth Doppler pass par's
 fp = 1                      ; first Doppler pass par's
 lp = 2                      ; last Doppler pass par's 

; default setup
if ~keyword_set(demo) then demo=0
if ~keyword_set(verbose) then verbose=0
if ~keyword_set(lm) then lm='mpf'
if ~keyword_set(avg) then avg=0
if ~keyword_set(cdnear_name) then cdnear_name=0

if keyword_set(nso) then iss_nm = 'nso' 
if keyword_set(nso) then iss_bc = 0.0d 

rootfile=starname+'_'+obsnm  ; used for naming the output files
savfile1='cd'+tag+'a'+rootfile   ;chunk_arr structure for 1st pass
savfile2='cd'+tag+'b'+rootfile   ;chunk_arr structure for 2nd pass
vdnm='vd'+tag+rootfile

c_light = dopenv.c_light
; COLORS FOR PLOTS
loadct,39, /silent
!p.background=255
!p.color=1

; GRAB VARIABLES FROM THE DOPENV STRUCTURE
; RENAME COMMON PATHS, FILES, INITIAL PSF DESCRIPTION
 files=dopenv.files_dir      ; e.g., /tous/mir1/files/

;old CCD mask:
;restore,dopenv.pfilt_file 

;create the ccd mask to down weight bad regions:
ccd_mask = dop_ccd_mask(dopenv=dopenv)

;now to mask out regions with telluric contamination:
ccd_mask = dop_telluric_mask(dopenv=dopenv, mask=mask)

if (pass eq 1) then begin
   print,'********************** FIRST PASS **************************'

  ; RETURN INITIAL CHUNK STRUCTURE 
  ; INITIAL WAVELENGTH SOLN AND PSF ARE FROM THE I2 OBSERVATION
  ; OPTIONAL: FITWAV WILL CALC CONTINUOUS WAVELENGTH SOLN FOR THE VD
  ; MATCH THE SYNTHETIC SPECTRUM TO OBSERVATION

	; kludge to handle old dsst's
	if dopenv.iss_nm ne 'iod' and dopenv.iss_nm ne 'nso' then begin  
	   if keyword_set(iss_nm) then restore, dopenv.files_dir+iss_nm ; restore the ISS
;          test=where(tag_names(dsst) eq 'W0',ntest)  & if ntest gt 0 then iss=dsst
	endif
	dop_chunk_setup, dopenv, chunk_arr, obs_info, $
		tag=tag, vdtag=vdtag, demo=demo, avg=avg, vdavgnm=vdavgnm, verbose=verbose, $
		crosscorl=crosscorl, tmpl_dir=tmpl_dir, yrmo=yrmo, cdnear_name=cdnear_name
 
    
	nchunks=n_elements(chunk_arr.ordt) 
  ; CH-LOOP: ONE CHUNK AT A TIME... 
   ;if keyword_set(yrmo) then avg=1 else avg=0
   for ch = 0, nchunks - 1 do begin
	  chunk = chunk_arr[ch]

	  ;chk_filt=filt[chunk.pixt:chunk.pixt+dopenv.n_pix-1, chunk.ordt]
	  chk_filt = ccd_mask[chunk.pixt:chunk.pixt+dopenv.n_pix-1, chunk.ordt]
	  ;stop
	  mod_chunk = dop_marq(chunk, dopenv, ch, pass=pass, chunk_arr=chunk_arr, $
						   demo=demo, verbose=verbose, lm=lm, smdisp=smdisp, $
						   smwav=smwav,pfilt=chk_filt, avg=avg) 
	; UPDATE CHUNK WITH MODEL PAR'S, VEL, FIT
	  chunk = mod_chunk   
	  chunk_arr[ch] = chunk
   endfor                    ; ch=CHUNK

   save,obs_info,chunk_arr,f=dopenv.files_dir+savfile1  
endif

if (pass eq 2) then begin 
   print,'*************** SECOND PASS - FREEZE PSF *********************'
   restore,dopenv.files_dir+savfile1
   nchunks=n_elements(chunk_arr.ordt)
   vd=dop_vdarr_setup(obsnm, dopenv, nchunks)

;      dop_fitwav,chunk_arr, npix=dopenv.n_pix, poly_ord=6, pass=2

   for i=0, nchunks - 1 do begin
	  chunk=chunk_arr[i]
	  chk_filt=filt[chunk.pixt:chunk.pixt+dopenv.n_pix-1, chunk.ordt]
	  mod_chunk = dop_marq(chunk, dopenv, i, pass=pass, chunk_arr=chunk_arr, $
		demo=demo, verbose=verbose, lm=lm, pfilt=chk_filt) 

	; UPDATE CHUNK AND VD_ARR WITH MODEL PAR'S, VEL, FIT
	  chunk = mod_chunk 
	  chunk_arr[i] = chunk
	  vd[i].ordt=chunk.ordt
	  vd[i].ordob=chunk.ordob
	  vd[i].pixt=chunk.pixt
	  vd[i].pixob=chunk.pixob
	  vd[i].w0=fix((*chunk.free_par)[lp].wcof[0])
	  vd[i].wcof[0]=(*chunk.free_par)[lp].wcof[0]-vd[i].w0
	  vd[i].wcof[1]=(*chunk.free_par)[lp].wcof[1]
	  vd[i].cts=chunk.cts
	  vd[i].weight=chunk.weight
	  vd[i].z=(*chunk.free_par)[lp].z
	  vd[i].vel=vd[i].z*c_light
	  vd[i].fit=chunk.fit[lp]
	  vd[i].psf=chunk.ip[*,lp]  
	  vd[i].par=(*chunk.free_par)[1].amp  ; zero for pass 2
	; END NLLS LOOP 
   endfor                
   save,obs_info, chunk_arr,f=dopenv.files_dir+savfile2
   save,vd, f=dopenv.files_dir+vdnm
   print,'Saving: ',savfile2
   print,'Saving: ',vdnm
endif 
chunk_clear, chunk_arr
;   winup,/all
end


