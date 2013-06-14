;+
;
;  NAME: 
;     dop_ccd_mask
;
;  PURPOSE: 
;   
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      dop_ccd_mask
;
;  INPUTS:
;
;  OPTIONAL INPUTS:
;
;  OUTPUTS:
;
;  OPTIONAL OUTPUTS:
;
;  KEYWORD PARAMETERS:
;    
;  EXAMPLE:
;      dop_ccd_mask
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2013.04.13 11:32:01
;
;-
function dop_ccd_mask, $
startpixel, $
endpixel, $
nord, $
dopenv = dopenv


rdsk, flat, dopenv.flatname
normflat = reform(flat[*,*,0])
flatdims = size(normflat, /dim)
ordlength = flatdims[0]

;now rotate 180 degrees to be the same orientation as the images:
normflat = rotate(normflat, 2)

mask = normflat * 0d + 1d
;chunk.pixt:chunk.pixt+dopenv.n_pix-1, chunk.ordt, dopenv=dopenv
;ord=0
;plot, normflat[*,ord], /xsty, /ysty, title=ord & ord++

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
end;dop_ccd_mask.pro