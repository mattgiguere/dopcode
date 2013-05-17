pro  psfav,vd,order,pixel,osamp,psf,nip,del_ord=del_ord,del_pix=del_pix, $
	  plot_key=plot_key,sigpsf=sigpsf,pixpsf=pixpsf,xra=xra,orddist=dummy

;common psfstuff,sigpsf,pixpsf,obpix
common psfstuff,psfsig,psfpix,obpix
;
;PURPOSE: Construct weighted average PSF within a domain on the echelle format
;INPUT: VD - structure
;       ORDER - order location at which PSF is to be determined
;       PIXEL - pixel location at which PSF is to be determined
;       OSAMP - 4 fine pixels within original
;OUTPUT:
;       PSF   - Final averaged PSF
;       NIP   - Number of IP's used in the domain to make the average
;
;KEYWORD: ORDDIST - positive integer - median pixel spacing between orders,
;                   this code assumes LICK standard 15 pixels per order,
;                   so the keyword ORDDIST is ignored
;
        if n_params() lt 1 then begin
          print,'Syntax:'
          print,'psfav,vd,ordt,pixt,osamp,psf,nip,del_ord=del_ord,' 
          print,'     del_pix=del_pix, plot_key=plot_key, '
;          print,'     sigpsf=sigpsf, pixpsf=pixpsf'
;          print,'i.e.,:'
;          print,'pixpsf=[0.00,-2.00,-1.25,-0.75, 0.75, 1.25, 2.00, 0.00, 0.00,0.00,0.00]'
;          print,'sigpsf=[0.00, 0.50, 0.40, 0.35, 0.35, 0.40, 0.50, 0.00, 0.00,0.00,0.00]'
          return
        end

tagnam=tag_names(vd)                                           ;VD tag_names
ordindex=first_el(where(tagnam eq 'ORDER'))                    ;order_index
if ordindex eq -1 then ordindex = first_el(where(tagnam eq 'ORDT'))   ;ORDER or ORDT?

	if n_elements(xra) ne 2 then xra=[-4,4]                ;plotting region
	if n_elements(sigpsf) gt 0 $
	   and n_elements(sigpsf) eq n_elements(pixpsf) then begin
;	     print,'Using a non-standard PSF in PSFAV'
	     psfsig=sigpsf
	     psfpix=pixpsf
        endif
        if n_elements(plot_key) ne 1 then plot_key=0
	xpsf=findgen(30*osamp+1)/osamp-15.
	if n_elements(del_pix) ne 1 then del_pix=100
	if n_elements(del_ord) ne 1 then begin
	   del_ord=4
	   if order le min(vd.(ordindex)) then del_ord=5
	   if order ge max(vd.(ordindex)) then del_ord=5
        endif
	ipind = where(abs(vd.(ordindex)-order) le del_ord and $
		    abs(vd.pixt-pixel) le del_pix)
        psfpar = vd(ipind)    ;subset of vd, boxed by del_ord and del_pix

; Take S/N into account,   PB 4/14/95
	psffit = psfpar.fit1 * sqrt(psfpar.cts)
	if n_elements(psfpar) eq 1 then ipfit=psffit $
	   else ipfit =  median(psffit)
	psfpar = psfpar(where(psffit lt (2.*ipfit)))  ;select those w/i 2*median psffit
;  PB moved the following two lines down 4 lines to fix bug, Sept. 13, 1994
	psflot=fltarr(n_elements(xpsf),n_elements(psfpar))
;	psfpar = psfpar(where(psfpar.pixt ne 540)) kludge, don't use pixel 540
;   Vertical Dist.; ~15 pxl/order
	ordist = (abs(order-psfpar.(ordindex)))*15
;   Horiz distance in pixels
	pxdist = fix( abs(pixel-psfpar.pixt) )
        ordwt = 1. / sqrt(ordist/40.+1)
        pxwt =  1. / sqrt(pxdist/40.+1)
;   WEIGHT for IP averaging
	;psfwt = ( (1./psffit)^2. ) * ordwt * pxwt
;GM -->  psfwt made proportional to 1/fit, not 1/fit^2
	psfwt = ( 1./psffit ) * ordwt * pxwt
;                   fit          order weight pixel weight  
        psfwt=psfwt/total(psfwt)
        inpr=xpsf*0.                       ;zero the psf
        numpsf = n_elements(psfwt)
	FOR qq = 0,numpsf-1 do begin

	     dumpsf = gpfunc(xpsf,float(psfpar(qq).par))
;GM -->    Fit parabola to 5 pixels around x = 0. --- Force centering.
             hm = 0.5 * max(dumpsf)                            ;half max
             ind  = where(dumpsf ge hm and abs(xpsf) lt 2.,np) ;Use Peak
             IF np ge 3 then begin             
               coef = pl_fit( xpsf(ind), dumpsf(ind), 2)
               cent = -0.5*coef(1)/coef(2) ;& print,cent       ;PSF Center
;PB Kludge, toss horrible PSF's and prevent BOMBS,  Dec. 6, 1993
;	       if abs(cent) lt 0.6 then begin ;Shift PSF
;PB Kludge Kludge, December 28, 1995   HR 2047
	       if abs(cent) lt 1.2 then begin ;Shift PSF
	       dumpsf = gpfunc(xpsf + cent, float(psfpar(qq).par))
	       endif else begin & psfwt(qq)=0. & dumpsf=xpsf*0. & endelse
             ENDIF

;GM -->    End Centering of PSF
;PB Kludge, toss horrible PSF's and prevent BOMBS,  Nov. 18, 1995
           if n_elements(dumpsf) ne n_elements(xpsf) then begin
	      psfwt(qq)=0.  &  dumpsf=xpsf*0.
           endif
	   psflot(*,qq) = dumpsf
	   inpr = inpr + psfwt(qq)*dumpsf
        ENDFOR
	inpr = osamp*inpr/total(inpr)  ;Normalize
	psfdif = psfwt  * 0.
	diflot = psflot * 0.
;plot,xpsf,inpr,xr=[-4,4]

	FOR qq = 0,numpsf-1 do begin
	   ;diflot(*,qq) = diflot(*,qq)-inpr
	   diflot(*,qq) = psflot(*,qq) - inpr
	   ;psfdif(qq) = total(abs(diflot(*,qq)))
	   psfdif(qq) = max((diflot(*,qq)))
        ENDFOR

	if n_elements(psfdif) eq 1 then medif = psfdif $
;	   else medif = median(psfdif)
	   else medif = mean(psfdif)
	spam = where(psfdif gt (1.1*medif),nspam)
	if nspam gt 0 then psfwt(spam)=0.        ;No weight to poor PSF's
	psfwt = psfwt/total(psfwt)
        inpr  = xpsf*0.                          ;zero the psf

 numgd = n_elements(where(psfwt gt 0.))  ;number of good psfs    gm,ew
 ;Averaging Section  (GM 1993,Oct.2 --- killed by GM 1995, Nov.12)
;        IF numgd gt 3 then BEGIN
;          Linear Fit: PSF(i) vs. pixel  (weighted by psfwt)
;           x = psfpar.pixt                             ;horiz posn of PSFs
;           FOR i = 0,n_elements(xpsf)-1 do begin       ;loop thru pixels
;             coef = polyfitw(x,psflot(i,*),psfwt(*),1) ;fit line to PSF 
;             inpr(i) = inpr(i) + poly([pixel],coef)    ;Eval. at pixel
;           END
;        END ELSE BEGIN
;          Median PSF Profiles (oldstyle)
psfind = where(psfwt gt 0)
           FOR i = 0,n_elements(xpsf)-1 do begin    ;loop thru pixels gm,ew
;The next line is a critical BUG, JESUS MOTHERFUCKING CHRIST! Dec 23, 1995
;               inpr(i) = median(psflot(i,*))   ;take median  gm,ew
               inpr(i) = median(psflot(i,psfind))
           END ;end median determination
;        END

;End Averaging 

        inpr = osamp*inpr/total(inpr)
	nip = n_elements(where(psfwt gt 0))

	IF plot_key eq 1 then begin
           !p.thick=1
	   xt = '!5Pixel = '+strtrim(pixel,2)
	   tt = '!5Order = '+strtrim(order,2)
           xcen = osamp * 15
           xsc = findgen(osamp*30)/osamp - 15.
	   yr=[0.,max(inpr)+.1]
	   plot,xsc,inpr,xr=xra,/xsty,symsize=.5,yra=yr,title=tt,xtitle=xt,charsize=1.8
	   FOR qq=0,(n_elements(psfwt)-1) do begin
	     if psfwt(qq) eq 0 then begin
                oplot,xsc,psflot(*,qq),co=191,ps=7,symsize=.1    ;bad
             end else begin
               oplot,xsc,psflot(*,qq),co=70+qq*20,ps=7,symsize=.6 ;good
             end
           ENDFOR ; qq
	   oplot,xsc,inpr,thick=3
        ENDIF  ;plot_key

	psf=inpr

return

end

