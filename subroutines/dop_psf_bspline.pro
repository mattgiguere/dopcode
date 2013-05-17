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

;the inverse variance. FIX THIS:
;invvarr = 1d / sqrt(10d * (psfarr + 1d))

;create the breakpoint array:
;places = [0,30,45,50,53,56,59,61,64,67,70,75,90,120]
places = dopenv.psfbsplnplaces
fullbkpt = bspline_bkpts(xarr, nord=bsord, placed=places)

;create the structure that will be used to generate the BSPSF:
sset = create_bsplineset(fullbkpt, bsord)

;now set the bspline coefficients to the values given from MPFIT.
;The range cuts off any additional, unnecessary parameters passed in:
sset.coeff = psfpar[0:n_elements(sset.coeff)-1]

;finally generate the PSF model that will be returned:
psfmodel = bspline_valu(xarr, sset)

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

