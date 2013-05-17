function lick_init, obsnm, iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag

; fischer oct 2009
; dopenv=lick_init('ri37.53','iod',0.0,iss_obnm='ri37.53')

; SETUP YOUR DIRECTORIES AND RENAME THE PATHS BELOW
  root_path='/home/fischer/dop2/'            ; the code lives here
  obs_dir = '/tous/mir1/fitspec/'                 ; the fits files (obs) live here

 ; for later - the initial guesses are stored in the "ipcf" file
  ipcf_dir = root_path+'ipcf_files/' ; the ipcf directory
  ipcf_st = 'ipcf_lick'+tag+'.dat' 
;  ipcf_st = 'ipcf_lickd.dat' 
  ipcf_file = ipcf_dir+ipcf_st

 ; the FTS spectrum is here
  fts_dir = root_path+'atlas/'
;  fts_atlas = 'ftseso50.bin'    ; fts_lick05.bin gives poor fits??
  fts_atlas = 'ftslick05.bin'    
;  fts_atlas = 'iodine_lickq1_pnnl.sav' ; very high sampling
;  fts_atlas = 'lickq1_pnnl_red10.sav'   ; sampling reduced by factor of ten

 ; the barycentric correction file is here 
  bary_file = '/tous/mir1/bary/bcvel.dat/'
  files_dir = '/tous/mir1/files/'    ; store the output (vd) files here

 ; 1851 pixels per order at Lick, start pixel = 50, end pixel = 1809
 ; iodine starts in order 54, at 5020 A, pixel 1090 to 1809 (9 chunks)
 ; good iodine lines in orders 35 - 53 (22 chunks per order = 418 chunks)
 ; useable iodine in order 34 from pixel 50 to 1170 = 14 chunks 
 ; total number of chunks: 441 

 ; 2048 pixels per order at CTIO, start pixel=520, end pixel=1479 
 ; iodine starts in order 5013 A, 12 chunks per order 
 ; useable orders 5 - 23, 12 chunks/ord * 18 orders = 216 chunks

  st_pix = 50   ; 
  st_ord = 38   ; 
  n_pix = 80
  n_ord = 16     
;  nch_ord = 22   pre-dewar4
  nch_ord = 27 
  n_chunks= n_ord*nch_ord    
  n_wcof = 2                 ; number of coefficients for fitting wavelength of one chunk

  ; I store Lick stuff in /mir1, hence variable names for the PSF parameters
  ;           [ 1      2   3     4     5     6      7     8    9    10    11   12    13    14   15   16  17]
  m1_ampl =  [0.25,  0.0,0.001,0.003,0.006, 0.03, 0.11, 0.21, 0.21, 0.11,0.03,0.007,0.002,0.001,0.0, 0.0, 0.0]
; old 
;  m1_psfpix = [0.0, -6.5, -5.5, -4.5, -3.5, -2.5, -1.8, -1.2,-0.6, 0.6, 1.2,  1.8,  2.5,  3.5, 4.5, 5.5, 6.5]
;  m1_psfsig = [1.4,  1.2,  1.2,  1.2,  1.2,  0.9,  0.8,  0.8, 0.7, 0.7, 0.8,  0.8,  0.9,  1.2, 1.2, 1.2, 1.2]
;  m1_psfpix = [0.0,-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0,-0.5, 0.5, 1.0,  1.5,  2.0,  2.5, 3.0, 3.5, 4.0]
;  m1_psfsig = [0.9, 1.2,  1.0,  0.9,  0.8,  0.8,  0.8,  0.7, 0.6, 0.6, 0.7,  0.8,  0.8,  0.8, 0.9, 1.0, 1.2]
  m1_psfpix = 1.1*[0.0,-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0,-0.5, 0.5, 1.0,  1.5,  2.0,  2.5, 3.0, 3.5, 4.0]
  m1_psfsig = [1.0, 1.2,  1.0,  0.9,  0.8,  0.8,  0.8,  0.7, 0.6, 0.6, 0.7,  0.8,  0.8,  0.8, 0.9, 1.0, 1.2]
  psfpix = m1_psfpix
  psfsig = m1_psfsig
  psf_ampl = m1_ampl 
  osamp = 4

 ; del_ord and del_pix are used for smoothing the PSF in the second
 ; pass - this suggests that the PSF is constant over this spatial 
 ; size on the CCD.  Here, I'm guessing it is smooth over 3
 ; orders (cross-disp) and 120 pixels (disp direction)
  del_ord = 3
  del_pix = 160                
  ord_dist = 15  ;typical separation between orders in pixels

; for SBs it can be useful to cross correlate instead of using 
; the barycentric shift
  xcorl_lambda = 4976.
  xcorl_pix1 = 1000 
  xcorl_npix = 120  
  xcorl_ord = 55    ;[findgen(500)+900,3]
  xcorl_shiftrange=18     ; make this larger for SB's - otherwise used to refine the wav guess
  decker = 'fiber'  ; "B5" or "B1" for Keck, "64" or "80" for Lick
  file_ext = 0

;stop
  dopenv = dop_init(obsnm, iss_nm, iss_bc, iss_obnm = iss_obnm, bary_file = bary_file, del_ord = del_ord, $
                del_pix = del_pix, ord_dist = ord_dist, npix = n_pix, files_dir = files_dir, obs_dir=obs_dir, $
                fts_dir = fts_dir, fts_atlas = fts_atlas, ipcf_file = ipcf_file,$
                psfpix = psfpix, psfsig = psfsig, psf_ampl = psf_ampl, osamp=osamp, $
                st_pix=st_pix, st_ord=st_ord, n_ord=n_ord, nch_ord = nch_ord, n_chunks=n_chunks, $ 
                decker = decker, dpad = dpad, n_wcof = n_wcof, xcorl_lambda = xcorl_lambda, $
                xcorl_pix1=xcorl_pix1, xcorl_ord=xcorl_ord, xcorl_shiftrange=xcorl_shiftrange,xcorl_npix=xcorl_npix,$
                file_ext = file_ext)


  return,dopenv
end
