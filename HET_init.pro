function het_init, obsnm, iss_nm, iss_bc, iss_obnm=iss_obnm

   root_path = '/Users/jtwright/idl/dop/'
   obs_dir = '/Users/jtwright/Data/HET/reduced/'

   ipcf_dir = root_path+''
   ipcf_st = 'ipcf.dat' 
   ipcf_file = ipcf_dir+ipcf_st


;   fts_dir = '/mir1/atlas/'
   fts_dir = '/Users/jtwright/Data/iodine/'
   fts_atlas = 'ftseso50.bin'    ;Scan of iodine cell

;   bary_file = '/mir1/bary/hbcvel.dat/'
   bary_file = '/Users/jtwright/Data/HET/hbcvel.dat'
   files_dir = root_path+'vd/'


   st_pix = 1000  ;first pix
   st_ord = 2   ;fird ord
   del_ord = 2    ; for PSF smoothing
   del_pix = 175                ; for PSF smoothing
   ord_dist = 65     ; physical distance for PSF smoothing  -- not used?  fix
   n_pix = 80        ;
   nch_ord = 38
   n_ord = 18
   nchunks = (n_ord-2)*nch_ord+30  ;orders 2 and 19 are short
   n_wcof = 2   ; Order of wavelength solution in a chunk

   file_ext = 2

   xcorl_lambda = 4976
;   m1_psfpix = 1.2*[0.00, -2.2,-1.7, -1.2, -0.7,-0.3, 0.3, 0.7, 1.2, 1.7, 2.2]
;   m1_psfsig = 1.2*[0.60,  0.6, 0.6,  0.6,  0.6, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6]
;   m1_ampl = m1_psfsig*0+0.02
  
   het_psfpix = [0.0, -4.5, -3.8, -3.2, -2.6, -2.0, -1.5, -0.9, -0.5, 0.5, 0.9, 1.5, 2.0, 2.6, 3.2, 3.8, 4.5]
   het_psfsig = [1.85, 0.8,  0.8,  0.7,  0.7,  0.6,  0.6,  0.5,  0.4, 0.4, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8]
   het_ampl =  [0.25,  0.0, 0.001, 0.003, 0.006, 0.03, 0.11, 0.21, 0.21, 0.11, 0.03, 0.007, 0.002, 0.001, 0.0, 0.0, 0.0]
   osamp = 4

   psfpix = het_psfpix
   psfsig = het_psfsig
   psf_ampl = het_ampl

;   psfpix = [0.00, -3.50, -2.90, -2.40, -1.90, -1.50, -1.20, -0.90, 0.90, 1.20, 1.50, 1.90, 2.40, 2.90, 3.50]
;   psfsig = [1.00, 0.00, 1.00, 0.80, 0.65, 0.50, 0.40, 0.30, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00, 0.00]
   decker = '60k'
   xcorl_limbda = 4967  ; find wavelength solution here later.

  dopenv = dop_init(obsnm, iss_nm, iss_bc, iss_obnm = iss_obnm, bary_file = bary_file, del_ord = del_ord, $
                del_pix = del_pix, ord_dist = ord_dist, npix = n_pix, files_dir = files_dir, obs_dir=obs_dir, $
                fts_dir = fts_dir, fts_atlas = fts_atlas, ipcf_file = ipcf_file,$
                psfpix = psfpix, psfsig = psfsig, psf_ampl = psf_ampl, osamp=osamp, $
                st_pix=st_pix, st_ord=st_ord, n_ord=n_ord, nch_ord = nch_ord, n_chunks=n_chunks, $ 
                decker = decker, dpad = dpad, n_wcof = n_wcof, xcorl_lambda = xcorl_lambda, $
                file_ext = file_ext)

  return, dopenv

end
