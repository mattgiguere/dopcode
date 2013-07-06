function ctio4k_init, obsnm, iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag, $
	date=date, dpad=dpad, avg=avg, psfmod=psfmod, shft_style=shft_style

; fischer oct 2009
; DF revise 13 Mar 2011 for e2v-4k ccd 
; dopenv=ctio_init('rqa10.7353','iod',0.0,iss_obnm='rqa10.7353')

date = strt(date)
year = '20'+strmid(date, 0, 2)

; SETUP YOUR DIRECTORIES AND RENAME THE PATHS BELOW
;  root_path='/tous/mir7/dopcode/dop/'            ; the code lives here
root_path='/tous/mir7/dop/'
if file_test('/Users/matt/') then root_path='/Users/matt/projects/dopcode/'
if file_test('/Users/debrafischer/') then root_path='/Users/debrafischer/dop1/'
if date eq '11' then obs_dir = '/tous/mir7/fitspec/' $
else obs_dir = '/tous/mir7/fitspec/'+date+'/'                         ; the fits files (obs) live here
observatory='ctio4k'

; for later - the initial guesses are stored in the "ipcf" file
ipcf_dir = root_path+'ipcf_files/' ; the ipcf directory
;;fischer jan7, 2013  ipcf_st = 'ipcf_ctio'+tag+'.dat' 
ipcf_st = 'ipcf_ctio'+tag+'.dat'
ipcf_file = ipcf_dir+ipcf_st

; the FTS spectrum is here
fts_dir = root_path+'atlas/'
fts_atlas = 'ctiov1f_pnnl_red5.dat'  ; V1 cell at CTIO
;  if strlen(tag) eq 2 then if strmid(tag,1,1) eq 's' then 
;  fts_atlas = 'ftsas4.bin'              ; McMath FTS for V1 cell 
;   fts_atlas = 'iodine_ctiov1_pnnl.dat' ; too highly sampled V1 cell at CTIO
;  if strlen(tag) eq 1 then fts_atlas = 'ftseso50.bin'  ; dummy fts (old lick) 

; pixel filter 
pfilt_file = 'torrent_filt.dat'

;pixel filter superceded by dop_ccd_mask.pro 130429 ~MJG
flatname = '/tous/mir7/flats/chi'+date+'.narrow.flat'

;add telluric mask 130614 ~MJG
;first restore the logstructure to search for one:
telids = -1
;either specify a filename, or loop through until you get a non-I2 B star:
telluricfn = ''
telluricfn = '/tous/mir7/fitspec/130506/achi130506.1180.fits'
loopdate = date
jddate = julday(strmid(loopdate, 2,2), strmid(loopdate,4,2), long('20'+strmid(loopdate, 0,2)))
while telids lt 0  and telluricfn eq '' do begin
  restore, '/tous/mir7/logstructs/'+year+'/'+loopdate+'log.dat'
  ;then the telluric file indeces are:
  telids = where(strt(log.decker) eq 'narrow_slit' and strt(log.propid) eq '80' and strt(log.iodcell) eq 'OUT')
  if telids lt 0 then begin
	;try one day earlier if need be:
	jddate--
	loopdate = jul2cal(jddate, /yymmdd)
  endif
endwhile
if telluricfn eq '' then $
	telluricfn = '/tous/mir7/fitspec/'+loopdate+'/achi'+$
	loopdate+'.'+log[telids].seqnum+'.fits'

; the barycentric correction file is here 
bary_file = '/tous/mir7/bary/qbcvel.dat'
files_dir = '/tous/mir7/files_df/'    ; store the output (vd) files here

; 2048 pixels per order at CTIO, start pixel=520, end pixel=1479 
; iodine starts in order 5013 A, 12 chunks per order 
; useable orders 5 - 23, 12 chunks/ord * 18 orders = 216 chunks

;TAG DEFINITIONS:
;0th element: chunk size & reduction code used
;    'd': 160 pixel chunks and old reduction
;    'a': 80 pixel chunks and Andrei reduction and Bstar I2s
;    'b': 80 pixel chunks, Andrei reduction, and Quartz-I2 only
;    'k': 80 pixel chunks and old reduction
;1st element: decker position
;    'a,c,d,f': narrow slit ("a" is COM PSF, "b" is bspline, "c" is peak-centering PSF)
;    'r': regular slit
;    's': slicer
;3rd elemnt: 
;    'f' smaller chunk averaging (17 weighted chunks)
;    'g' even smaller chunk averaging (9 weighted chunks)

st_pix = 80
;  st_ord = 12    ; could pick up half of order 11...
st_ord = 13    ; Torrent controller Mar 2012
npsfpix = 121

if strmid(tag,0,1) eq 'd' then begin
  n_pix = 160  ;18
  nch_ord = 19 ;38             ; at CTIO, number of chunks per order 
endif
if strmid(tag,0,1) eq 'a' or strmid(tag,0,1) eq 'b' then begin
 ;NOTE: ALSO NEED TO EDIT THE BOTTOM OF THIS FILE FOR ANDREI'S REDUCTION TO WORK:
  n_pix = 80  ;18
  nch_ord = 38 ;38             ; at CTIO, number of chunks per order 
endif
if strmid(tag,0,1) eq 'k' then begin
  n_pix = 80  
  nch_ord = 38                 ; at CTIO, number of chunks per order 
endif
if strmid(tag,0,1) eq 'n' then begin
;20120419
  n_pix = 80  
  nch_ord = 38                 ; at CTIO, number of chunks per order 
endif
if strmid(tag,0,1) eq 't' then begin
  n_pix = 80  
  nch_ord = 38                 ; at CTIO, number of chunks per order 
endif
if strmid(tag,0,1) eq 's' then begin
  n_pix = 80  
  nch_ord = 38                 ; at CTIO, number of chunks per order 
endif
if tag eq 'rchi' then begin
  n_pix = 80  
  nch_ord = 38                 ; at CTIO, number of chunks per order 
endif
n_wcof = 2              ; number of coefficients for fitting wavelength of one chunk ([0, 1, 2])
n_ord = 20     
;  n_ord = 16  ; telluric contamination is significant in orders 29 and 30.
n_chunks= n_ord*nch_ord    ; 710 for keck

; not sure I like the following (DF 6/14/12)
;  decker = 'narrow'  ; CTIO: narrow, slit, slicer, fiber; Keck:B5 or B1; Lick: 64 or 80
;  if decker eq 'narrow' then scale = 0.8 else scale = 1.

; I store CTIO stuff in /mir7, hence variable names for the PSF parameters
; NOTE: in pass 0 of the doppler code, the FWHM of the central gaussian is
; a free parameter.  Therefore, m7_psfsig will be reset and then fixed.
;           [ 1    2    3     4     5     6     7     8    9     10    11   12    13   14    15   16   17]
m7_ampl =  [0.6,  0.0, 0.0,0.001,0.003,0.006, 0.03, 0.11, 0.21, 0.21, 0.11,0.03,0.007,0.002,0.001,0.0, 0.0]
m7_psfpix = [0.0,-4.0, -3.6, -3.0, -2.4, -1.8, -1.2, -0.8, -0.4, 0.4, 0.8,  1.2,  1.8,  2.4, 3.0, 3.6, 4.0]
m7_psfsig = [0.9, 0.0,  0.8,  0.7,  0.7,  0.7,  0.7,  0.6,  0.5, 0.5, 0.6,  0.7,  0.7,  0.7, 0.7, 0.8, 0.0]

psfpix = m7_psfpix
psfsig = m7_psfsig
psf_ampl = m7_ampl 
osamp = 4
if ~keyword_set(psfmod) then psfmod='gaussian'

; del_ord and del_pix are used for smoothing the PSF in the second
; pass - this suggests that the PSF is constant over this spatial 
; size on the CCD.  Here, I'm guessing it is smooth over 3
; orders (cross-disp) and 120 pixels (disp direction)
del_ord = 1  ;2;3
del_pix = 120 ;160;240                
ord_dist = 30  ;typical separation between orders in pixels

; for SBs it can be useful to cross correlate instead of using 
; the barycentric shift
xcorl_lambda = 4976.
xcorl_pix1 = 1300
xcorl_npix = 500
xcorl_ord = 10     ;[findgen(500)+900,3]
xcorl_shiftrange=20     ; make this larger for SB's - otherwise used to refine the wav guess
; decker = 'narrow'  ; "B5" or "B1" for Keck, "64" or "80" for Lick
file_ext = 0

dopenv = dop_init(obsnm, iss_nm, iss_bc, iss_obnm = iss_obnm, bary_file = bary_file, del_ord = del_ord, $
			  del_pix = del_pix, ord_dist = ord_dist, npix = n_pix, files_dir = files_dir, obs_dir=obs_dir, $
			  fts_dir = fts_dir, fts_atlas = fts_atlas, ipcf_file = ipcf_file, pfilt_file=pfilt_file, $
			  psfpix = psfpix, psfsig = psfsig, psf_ampl = psf_ampl, osamp=osamp, avg=avg, observatory=observatory, $
			  st_pix=st_pix, st_ord=st_ord, n_ord=n_ord, nch_ord = nch_ord, n_chunks=n_chunks, $ 
			  psfmod=psfmod, npsfpix=npsfpix, shft_style=shft_style, $
			  decker = decker, dpad = dpad, n_wcof = n_wcof, xcorl_lambda = xcorl_lambda, $
			  xcorl_pix1=xcorl_pix1, xcorl_ord=xcorl_ord, xcorl_shiftrange=xcorl_shiftrange,xcorl_npix=xcorl_npix,$
			  file_ext = file_ext, flatname = flatname, telluricfn=telluricfn)

if strmid(tag,0,1) eq 'a' or strmid(tag,0,1) eq 'b' then begin
 dopenv.obsnm = 'a'+strmid(dopenv.obsnm, 1, strlen(dopenv.obsnm)-1)
endif
if strmid(date,0,2) eq '11' then begin
  dopenv.obsnm = 'r'+strmid(dopenv.obsnm,1,strlen(dopenv.obsnm)-1)
endif 

return,dopenv
end
