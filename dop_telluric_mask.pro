;+
;
;  NAME: 
;     dop_telluric_mask
;
;  PURPOSE: 
;   
;
;  CATEGORY:
;      DOPCODE
;
;  CALLING SEQUENCE:
;
;      ccd_mask = dop_telluric_mask(dopenv=dopenv)
;
;  INPUTS:
;		DOPENV: The Doppler Code environmental structure
;
;  OPTIONAL INPUTS:
;		MASK: Input this if there is already an array of the proper 
;			dimensions that you want to start with, otherwise it will
;			start with a clean slate. This was intended for input from
;			the mask that was already created from DOP_CCD_MASK that 
;			masks out bad regions on the detector.
;
;  OUTPUTS:
;
;  OPTIONAL OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      ccd_mask = dop_telluric_mask(dopenv=dopenv, mask=mask)
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2013.06.13 12:06:57
;
;-
function dop_telluric_mask, $
dopenv=dopenv, $
mask=mask

angstrom = '!6!sA!r!u!9 %!6!n'
!p.color=0
!p.background=255
loadct, 39, /silent
usersymbol, 'circle', /fill, size_of_sym = 0.5

;read in the B star spectrum that'll be used for the 
;telluric mask:
telim = readfits(dopenv.telluricfn, hd)
telspec = reform(telim[1,*,*])
telwav = reform(telim[0,*,*])

;now read in the flat:
rdsk, flat, dopenv.flatname
normflat = reform(flat[*,*,0])
flatdims = size(normflat, /dim)
ordlength = flatdims[0]

;now rotate 180 degrees to be the same orientation as the images:
normflat = rotate(normflat, 2)
flatmod = rotate(reform(flat[*,*,1]),2)

stop
if ~keyword_set(mask) then mask = normflat * 0d + 1d

;orders 0-*: 
threshold=0.5
for ord=0, n_elements(mask[0,*])-1 do begin
  x = where(normflat[*,ord] le threshold, xct)
  ;now to go through and pad the sides of bad areas:
  for i=0, xct-1 do begin
  ;the left and right padding amount:
  lpad = 3 & rpad = 3
  ;special padding for areas close to edge:
  if x[i] lt lpad then lpad = x[i]
  if x[i] gt (ordlength - rpad - 1) then rpad = ordlength - x[i] - 1
  mask[(x[i]-lpad):(x[i]+rpad),ord] = 0d
  endfor
  
endfor

;restore, dopenv.pfilt_file
;for ord=0, 61 do begin
;  usersymbol, 'circle', size_of_sym=1, /fill
;  plot, filt[*,ord], ps=8, /xsty, yran=[-0.1,1.1], title=ord
;  usersymbol, 'circle', size_of_sym=0.5, /fill
;  oplot, (normflat[*,ord]-0.9)*10d, col=75, ps=8
;  oplot, mask[*,ord], col=250, ps=8
  ;stop
;  endfor
;stop
return, mask





stop
end;dop_telluric_mask.pro