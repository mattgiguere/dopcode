;+
;
;  NAME: 
;     dop_psf_bspline
;
;  PURPOSE: This function creates a PSF model using bsplines
;   
;
;  CATEGORY:
;      CHIRON
;
;  CALLING SEQUENCE:
;
;      dop_psf_bspline
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
;      dop_psf_bspline
;
;  MODIFICATION HISTORY:
;        c. Matt Giguere 2013.03.13 18:03:31
;
;-
function dop_psf_bspline, $
	psfpar, $
	dopenv=dopenv, $
	xarr=xarr, $
	cntr=cntr

if ~keyword_set(xarr) then xarr = dindgen(121)

bsord = dopenv.psfbsord

;create the breakpoint array:
places = dopenv.psfbsplnplaces
fullbkpt = bspline_bkpts(xarr, nord=bsord, placed=places)

;create the structure that will be used to generate the BSPSF:
sset = create_bsplineset(fullbkpt, bsord)

;now set the bspline coefficients to the values given from MPFIT.
;The range cuts off any additional, unnecessary parameters passed in:
sset.coeff = psfpar[0:n_elements(sset.coeff)-1]

;finally generate the PSF model that will be returned:
psfmodel = bspline_valu(xarr, sset)

;Center the PSF model to overcome degeneracies with the 
;wavelength solution:
if keyword_set(cntr) then begin
	;half the # of elements in the PSF model:
	hpsf = n_elements(psfmodel)/2d

  ;the simple sort and median shift center:
  if cntr eq 2 then begin
	;sort the positions, making the position of the highest element
	;listed first and in descending order:
	arr = reverse(sort(psfmodel))

	;take the median position of the top ten percent:
	mmaxes = median(arr[0:ceil(0.05*n_elements(psfmodel))])

	;now shift by that median max amount towards the center:
	shftamount = hpsf - mmaxes
	psfmodel = shift(psfmodel, shftamount)
  endif ;cntr=2
  
  ;now the center of mass shift:
  if cntr eq 1 then begin
	totpsf = total(psfmodel)
	psfsum = dblarr(n_elements(psfmodel))
	psfsum[0]=psfmodel[0]
	;for i=1, n_elements(psfmodel)-1 do psfsum[i] = psfsum[i-1]+psfmodel[i]
	i=1L
	compos=0d
	while compos le 0 do begin
	  psfsum[i] = psfsum[i-1]+psfmodel[i]
	  ;store the position of the halfway point and exit:
	  if psfsum[i] ge totpsf/2d then compos=i
	  i++
	endwhile

	;now shift by the com towards the center:
	shftamount = hpsf - compos
	;loadct, 39
	;print, totpsf, psfsum[i-1], compos, i
	;plot, psfmodel
	psfmodel = shift(psfmodel, shftamount)
	;oplot, psfmodel, col=250
	;stop
  endif;cntr eq 1
endif;cntr
;if max(psfpar) gt 1 then stop


;if max(psfmodel) gt 1 then psfmodel = psfmodel/max(psfmodel)
;psfmodel = psfmodel/max(psfmodel)
psfmodel = psfmodel/total(psfmodel)
return, psfmodel
end;dop_psf_bspline.pro

;*******************************************************
; END OF CODE, JUST NOTES:
;*******************************************************
;Best fit to narrow slit PSF pars:
psfpar = [ $
  -0.00528339, $
   0.00141147, $
 -0.000838449, $
  0.000675785, $
   0.00221572, $
    0.0260462, $
     0.197992, $
     0.570373, $
     0.532363, $
    0.0704252, $
    0.0239822, $
   0.00139637, $
   0.00186412, $
  -0.00140414, $
   0.00205050, $
  -0.00716207 $
  ]

;for testing purposes:
if ~keyword_set(dopenv) then begin
psfpix = dindgen(121)
dopenv = {dopenv,                  $
	  psfmod: 'bspline', $;PSF model to use: either 'gaussian' or 'bspline'
	  psfbsplnplaces: [0,30,45,50,53,56,59,61,64,67,70,75,90,120], $ ;coefficients for the bspline model
	  psfbsinvvar: dblarr(n_elements(psfpix)), $ ;inverse variance for bspline weighting
	  psfbsord: 4 $ ; the order to use for the bspline
	 }
endif

