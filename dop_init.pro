;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_INIT
;
; Initial coding: David Abouav Aug 2006
; Revised: D. Fischer Jan 2008 
; bookkeeping revisions with J. Wright Oct 2009
;
; CALLED FROM:
;   DOP_MAIN
;
; PURPOSE: 
;   Establishes the Doppler code paths and other defaults 
;
; INPUTS:
;   TAG: observatories (and settings) are identified by the run tags
;   
;
; OUTPUTS:
;   DOPENV: structure with relevant paths and parameters
;
; OUTSTANDING ISSUE: add number of orders in wavelength soln
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_INIT, obsnm, iss_nm, iss_bc, iss_obnm = iss_obnm, $
				bary_file = bary_file, del_ord = del_ord, $
                del_pix = del_pix, ord_dist = ord_dist, npix = n_pix, $
                files_dir = files_dir, obs_dir=obs_dir, $
                fts_dir = fts_dir, fts_atlas = fts_atlas, ipcf_file = ipcf_file, $
                pfilt_file=pfilt_file, $
                psfpix = psfpix, psfsig = psfsig, psf_ampl = psf_ampl, $
                osamp=osamp, avg=avg, observatory=observatory, $
                st_pix=st_pix, st_ord=st_ord, n_ord=n_ord, nch_ord = nch_ord, $
                n_chunks=n_chunks, $ 
                psfmod=psfmod, npsfpix=npsfpix, shft_style=shft_style, $
                decker = decker, dpad = dpad, n_wcof = n_wcof, xcorl_lambda = xcorl_lambda, $
                xcorl_pix1=xcorl_pix1, xcorl_ord=xcorl_ord, $
                xcorl_shiftrange=xcorl_shiftrange,xcorl_npix=xcorl_npix,$
                file_ext = file_ext

   restore, bary_file             ; Restores bcat
   xbc = where(bcat.obsnm eq obsnm, nxbc) 
   if nxbc eq 0 then begin
	  xbc = where(strt(bcat.obsnm) eq 'a'+strt(strmid(obsnm, 1, strlen(obsnm)-1)), nxbc) 
	  if nxbc eq 0 then begin
	  xbc = where(strt(bcat.obsnm) eq strt(strmid(obsnm, 1, strlen(obsnm)-1)), nxbc) 
	  if nxbc eq 0 then begin
            xbc = where(strt(bcat.obsnm) eq obsnm,nxbc)
            if nxbc eq 0 then xbc=where('r'+strmid(obsnm,1,strlen(obsnm)-1) eq bcat.obsnm,nxbc)
            if nxbc eq 0 then stop,'probably not in qbcvel.ascii'  
	  endif
   endif
   endif
   obs_bc = bcat[xbc].bc
   obj_nm = bcat[xbc].objnm
   if obj_nm eq 'iodine' then obj_nm = 'iod' 

 ; Declaration of the 'dopenv' structure. 
   dopenv = {dopenv,                  $
             obsnm: obsnm,                    $ ; observation ID, input to dop_init
             obs_dir: obs_dir,                $ ; observation dir
             obs_bc: obs_bc,                  $ ; BC of observation
             obj_nm: obj_nm,                  $ ; starname
             decker: '?',                     $ ; slit width, e.g. 'B1' 
             iss_nm: iss_nm,                  $ ; intrinsic stellar spectrum, input to dop_init
             iss_obnm: iss_obnm,              $ ; observation name of IS
             iss_bc: iss_bc,                  $ ; BC of intrinsic spectrum
             osamp: 4,                        $ ; oversampling for model
             avg:0,                           $ ; use avg psf for 0th pass
			 observatory:' ',                 $ ; which observatory? 
             st_pix: 0,                       $ ; first pixel
             st_ord: 0,                       $ ; first chunk
             nch_ord: 0,                      $ ; number of chunks per order
             n_ord: 0,                        $ ; number of orders
             del_ord: 0,                      $ ; psf averaging: number of orders
             del_pix: 0,                      $ ; psf averaging: number of pixels
             ord_dist: 0,                     $ ; order separation 
             dpad: 12.0,                      $ ; padding on iss chunks  
             c_light: 2.99792458d8,           $ ; speed of light
             n_wcof: 2,                       $ ; order of polynomial fit in wavelength
             n_pix: 80,                       $ ; number of pixels per chunk
             n_chunks: 710,                   $ ; number of chunks in vd structure 
             npsfpix: 121,                    $ ; number of chunks in vd structure 
             n_iter: 3,                       $ ; number of iterations 
             shft_style: 0,                   $ ; 1 for peak-centering SLSF and 2 for COM centering
             xcorl_lambda: 0.0d,              $ ; wavelength for xcorl to find initial z
             xcorl_pix1:900,                  $ ; pixels for xcorl to find initial z
             xcorl_ord:3,                     $ ; order for xcorl to find initial z
             xcorl_npix:500,                  $ ; order for xcorl to find initial z
             xcorl_shiftrange:15,             $ ; pixel range for xcorl to find initial z
             file_ext: 0,                     $ ; file extension for mrdfits
             files_dir: '',                   $ ; Location where vd files should be saved.
             fts_dir: '',                     $ ; Name of FTS file
             fts_atlas: '',                   $ ; Name of FTS file
             ipcf_file: '',                   $ ; name of ipcf file
             pfilt_file: '',                  $ ; pixel filter file 
             bary_file: bary_file,            $ ; Full path and name of appropriate bcvel.dat file.
             psf_ampl: fltarr(n_elements(psfpix)),  $ ; amplitude of Gaussians for PSF description
             psfpix: fltarr(n_elements(psfpix)),  	$ ; location of Gaussians for PSF description
             psfsig: fltarr(n_elements(psfpix)),  	$ ; width of Gaussians for PSF description
             psfmod: '', 					  $ ;PSF model to use: either 'gaussian' or 'bspline'
             psfbsplnplaces: [0,30,45,50,53,56,59,61,64,67,70,75,90,120], $ ;pix positions for the bspline model
             psfbsinvvar: dblarr(n_elements(psfpix)), $ ;inverse variance for bspline weighting
             psfbsord: 4, $ ; the order to use for the bspline
             psfcntr: 1 $ ;1 if you want to center the PSF, 0 not to center
            }
   
  ; SUPERCEDE DEFAULT VALUES
   if n_elements(del_ord) gt 0 then dopenv.del_ord = del_ord
   if n_elements(del_pix) gt 0 then dopenv.del_pix = del_pix
   if n_elements(ord_dist) gt 0 then dopenv.ord_dist = ord_dist
   if n_elements(npix) gt 0 then dopenv.n_pix = npix
   if n_elements(n_chunks) gt 0 then dopenv.n_chunks = n_chunks
   if n_elements(files_dir) gt 0 then dopenv.files_dir = files_dir
   if n_elements(fts_atlas) gt 0 then dopenv.fts_atlas = fts_atlas
   if n_elements(fts_dir) gt 0 then dopenv.fts_dir = fts_dir
   if n_elements(ipcf_file) gt 0 then dopenv.ipcf_file = ipcf_file
   if n_elements(pfilt_file) gt 0 then dopenv.pfilt_file = pfilt_file
   if n_elements(psfpix) gt 0 then dopenv.psfpix = psfpix
   if n_elements(psfsig) gt 0 then dopenv.psfsig = psfsig
   if n_elements(psf_ampl) gt 0 then dopenv.psf_ampl = psf_ampl
   if n_elements(psfmod) gt 0 then dopenv.psfmod = psfmod
   if n_elements(st_pix) gt 0 then dopenv.st_pix = st_pix
   if n_elements(st_ord) gt 0 then dopenv.st_ord = st_ord
   if n_elements(n_ord) gt 0 then dopenv.n_ord = n_ord
   if n_elements(nch_ord) gt 0 then dopenv.nch_ord = nch_ord
   if n_elements(osamp) gt 0 then dopenv.osamp = osamp
   if n_elements(npsfpix) gt 0 then dopenv.npsfpix = npsfpix
   if n_elements(shft_style) gt 0 then dopenv.shft_style = shft_style
   if n_elements(avg) gt 0 then dopenv.avg = avg
   if n_elements(observatory) gt 0 then dopenv.observatory=observatory
   if n_elements(decker) gt 0 then dopenv.decker = decker
   if n_elements(dpad) gt 0 then dopenv.dpad = dpad
   if n_elements(n_wcof) gt 0 then dopenv.n_wcof = n_wcof
   if n_elements(xcorl_lambda) gt 0 then dopenv.xcorl_lambda = xcorl_lambda
   if n_elements(xcorl_pix1) gt 0 then dopenv.xcorl_pix1 = xcorl_pix1
   if n_elements(xcorl_ord) gt 0 then dopenv.xcorl_ord = xcorl_ord
   if n_elements(xcorl_npix) gt 0 then dopenv.xcorl_npix = xcorl_npix
   if n_elements(xcorl_shiftrange) gt 0 then dopenv.xcorl_shiftrange = xcorl_shiftrange
   if n_elements(obs_dir) gt 0 then dopenv.obs_dir =  obs_dir
   if n_elements(file_ext) eq 0 then dopenv.file_ext = file_ext

   return, dopenv  

END








