;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PRO DOP_DSST
;
; PROCEDURES CALLED: 
;
; PURPOSE: 
;   Generate a deconvolved stellar template
;   Based on algorithms described in: 
;      Valenti et al. 1995 PASP 107 966
;      Modified version of dsst.pro (Marcy & Butler, 1995)
;
; PROCEDURE:
;   1) extract the wavelength soln and psf from bookend iodine chunks
;   2) extract (coadded) template spectra "spec" for those same chunks
;   3) find a "function," dspec, such that dspec equals the 
;      template chunk convolved with psf
;
; CALLING SEQUENCE: 
;   dop_dsst, starnm='achi120606.1125, coadded_spec='hd128620_achi120606.dat', $
;	          obsnm='chi120606.1125', date='120606', observatory='ctio4k', $
;             iss_tag='ax2', iodnm='HR5267', $
;			  bookend_iod_obsnm=['1113','1128','1129','1130'], $
;	 		  plot=plot
;
; INPUTS: 
;   COADDED_SPEC: the coadded and raycleaned spectrum array 
;   OBSNM: observation number to find entry in qbcvel.dat 
;       For coadded templates, just use one of the middle observation id's.
;   DATE: needed for the file path
;   OBSERVATORY: ctio4k 
; 	FILT_FILE: bad pixel file "filt" used to identify bad pixels for the chisq fit
; 	VD_TAG: tag used to build the iodine vd structures 
;   ISS_TAG: for filename of output dsst structure: "dsst128620[tag]_obsnm.dat"
;   IODNM: object name for bookend I2 observation (e.g., 'iod' or 'HR5267')
;   BOOKEND_IOD_OBSNM: array of iodine obsnm to provide the psf for deconvolution
;
; OPTIONAL INPUT: 
;   J_BVAL: bvalue for Jansson deconvolution (0.95 - 1.05, typically) 
;   JANSSON: carry out Jansson deconvolution with parabolic smoothing
;	MAX_VD_CHI: the max chisq fit for "good" bookend iodines 
;	DSST_PAD: (default is npix/4) 
;   PLOT: 1: spec and dst plots
;         2: PSF plots and residuals to nonlinear fits of zero point) 
;		  3: residuals to nonlinear fits of dispersion
;  	SMOOTH: smooth the wavelength and dispersion with a polynomial fit across the order
;
; OUTPUTS
;   DSST: a structure with model parameters:
;         ORDT            INT       13
;         PIXT            INT       80
;         W0              DOUBLE    5030.1909
;         W1              DOUBLE    0.0047160774
;         WEIGHT          FLOAT     0.00325524
;         DST             FLOAT     Array[416]
;
; Debra Fischer, Yale June 2012
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO DOP_DSST, starnm=starnm, coadded_spec=coadded_spec, obsnm=obsnm, date=date, $
	          observatory=observatory, pfilt_file=pfilt_file, vdtag=vdtag, $
              iss_tag=iss_tag, j_bval=j_bval, jansson=jansson, vd_array=vd_array, $
              iodnm=iodnm, bookend_iod_obsnm=bookend_iod_obsnm, vdavg=vdavg, $
              max_vd_chi=max_vd_chi, dsst_pad=dsst_pad, plt=plt, $
              smth=smth, mov=mov, fun=fun, ngrid=ngrid, smooth_val=smooth_val

;IDL> tag='adg'
;IDL> vd_array1 = 'vd'+tag+'HR8425_achi120710.'+['1137','1138','1139']
;IDL> vd_array2 = 'vd'+tag+'HR472_achi120710.'+['1143','1144','1145']
;IDL> vd_array = [vd_array1, vd_array2]
;IDL> dop_dsst, starnm='10700',coadded_spec='hd10700_achi120710.dat', date='120710', $
;               obsnm='achi120710.1143', observatory='ctio4k', vdtag='adg', $
;               vd_array=vd_array,iss_tag='adg' ;,/mov
;IDL> dop_dsst, starnm='10700',coadded_spec='hd10700_achi120710.dat', date='120710',$
;			    obsnm='achi120710.1143',observatory='ctio4k',vdtag='adg', $
;				iss_tag='adg',/vdavg, /mov

c = 2.99792458d8      

;if ~keyword_set(starnm) then starnm = '128620'
;if ~keyword_set(starnm) then starnm = '10700'
;if ~keyword_set(coadded_spec) then coadded_spec='hd128620_achi120606.dat' 
;if ~keyword_set(coadded_spec) then coadded_spec='hd10700_achi120710.dat' 
;if ~keyword_set(obsnm) then obsnm = 'achi120606.1114'  ;obsnm for
;128620 template 
;if ~keyword_set(obsnm) then obsnm = 'achi120710.1143'  ;obsnm for 128620 template 
;if ~keyword_set(date) then date='120606'
;if ~keyword_set(date) then date='120710'
;if ~keyword_set(observatory) then observatory = 'ctio4k'
;if ~keyword_set(vdtag) then vdtag='am'
;if ~keyword_set(iss_tag) then iss_tag='am'
if ~keyword_set(j_bval) then j_bval=0.97
;if ~keyword_set(iodnm) then iodnm='HR472'
if ~keyword_set(max_vd_chi) then max_vd_chi=1.5
if ~keyword_set(plt) then plt=0
if ~keyword_set(mov) then mov=0
if ~keyword_set(fun) then fun=0
if ~keyword_set(ngrid) then ngrid=20
if ~keyword_set(smooth_val) then smooth_val=2
if ~keyword_set(dsst_pad) then dsst_pad=12 ;integer (pixel index)
bpad=30                                    ;integer (pix index)

if j_bval gt 1.0 then stop, 'REALLY??  this will introduce ringing!'

if observatory eq 'ctio4k' then restore,'/tous/mir7/bary/qbcvel.dat'
; if the obsnm is given as "achi.." then trim off the first character
if strmid(obsnm,0,4) eq 'achi' then obsnm=strmid(obsnm,1,strlen(obsnm)-1)

x=where(bcat.obsnm eq obsnm,nx)
if nx eq 0 then stop,'Observation not found in qbcvel.dat'
iss_nm=starnm             ;bcat[x].objnm
iss_bc=bcat[x].bc
iss_jd=bcat[x].jd

; get the dopenv structure
dopenv=ctio4k_init(obsnm, iss_nm, iss_bc, iss_obnm=obsnm, tag=iss_tag, $
   	   date=date, dpad=dsst_pad)

;dopcode_path='/Users/debrafischer/research/dop2/'
dopcode_path='/Users/debrafischer/dop2/'
osamp=dopenv.osamp    ; e.g. 4 (integer)
npix=dopenv.n_pix     ; e.g. 80  (integer)
path=dopenv.obs_dir   ; e.g /tous/mir7/fitspec/120606/
dfd=dopenv.files_dir
n_chunks=dopenv.n_chunks  ; integer
n_ord=dopenv.n_ord        ; integer
st_ord=dopenv.st_ord      ; integer
nch_ord=dopenv.nch_ord    ; integer
pfilt_file=dopenv.pfilt_file ; string

if ~keyword_set(vd_array) and ~keyword_set(vdavg) then begin
	vdprefix='vd'+vdtag+iodnm+'_'+'achi'+date+'.'
	;vdprefix='vd'+vdtag+iodnm+'_'+'rchi'+date+'.'
	if keyword_set(bookend_iod_obsnm) then vd_array=vdprefix+bookend_iod_obsnm $
	;  else vd_array=vdprefix+['1113','1114','1116','1128','1129', '1130']
	;else vd_array=vdprefix+['1113','1114']
	;else vd_array=vdprefix+['1143','1144','1145']
	else vd_array=vdprefix+['1251','1252','1289','1290','1291']
	print,'The bookend iodine observations are: '
	print, vd_array
	ans='n'
	read,'Is this set of bookend iodines correct (y/n)? ',ans
	if ans eq 'n' then stop
endif 

if ~keyword_set(vdavg) then begin
;restore the bookend iodine vd's 
	num_i2=1   ; counting the number of bookend iodines with good fits
	restore,dfd+vd_array[0]

	if median(vd.fit) gt max_vd_chi then stop,'This VD has a high chisq fit'
	vdx=vd
	for i=1,n_elements(vd_array)-1 do begin
    ;print,vd_array[i]
		restore,dfd+vd_array[i]
		if median(vd.fit) lt max_vd_chi then begin
			num_i2+=1
			vdx=[vdx,vd]
		endif else print,'This bookend iodine was rejected: ',vd_array[i]
	endfor 
	if num_i2 ne n_elements(vd_array) then stop, 'One or more of the VDs had a high chisq fit'

n_psf=n_elements(vdx[0].psf)      	; number of array elements in the psf (typically 121)
vd_avpsf=dblarr(n_chunks, n_psf)  	; weighted ave psf for each chunk
vd_avwav = dblarr(n_chunks)       	; weighted mean wavelength soln for each chunk 
vd_avdisp = dblarr(n_chunks)      	; weighted mean dispersion soln for each chunk 
err_wav = dblarr(n_chunks)	      	; stddev of mean wavelength soln
wt=dblarr(n_chunks,num_i2)        	; weights for each chunk of each observation 
w0=fltarr(nch_ord)                	; w0 across an order of the averaged bookend
wcof1=dblarr(nch_ord)             	; wcof1 across an order of the ave bookend I2
wcof2=dblarr(nch_ord)             	; dispersion across an order of the ave bookend
wav0=dblarr(n_ord,nch_ord,num_i2)  	; chunk wavelength for each bookend observation
disp=dblarr(n_ord,nch_ord,num_i2) 	; chunk dispersion for each bookend observation
xip=findgen(30*osamp+1)/osamp - 15  ; xarr for oversampled PSF

; get the weighted average psf
	for i=0,n_chunks-1 do begin
		for j=0,num_i2-1 do begin
			if j eq 0 then indx=[i+(j*n_chunks)] $
		          else indx=[indx, i+j*n_chunks] 
			wt[i,j]=(1./vdx[indx[j]].fit)^2
		endfor ; j (chunks for averaging psf)
;print,indx, wt[i,*]
		for k=0,n_psf-1 do vd_avpsf[i,k]=total(wt[i,*]*vdx[indx].psf[k])/total(wt[i,*])
		  if plt eq 2 then begin
			str1='!6 PSF for chunk: '+strcompress(string(i),/rem)
			for j=0, num_i2-1 do begin
				if j eq 0 then plot, findgen(n_psf), vdx[indx[j]].psf[*],col=0, $
					thick=2, xra=[38,82],/xsty, yra=[0,0.6],/ysty, xtitl=str1
				if j gt 0 then oplot, findgen(n_psf), vdx[indx[j]].psf[*],col=j*45, $
					thick=2
			endfor
			!p.color=222
			xyouts,40,0.44,/data,'!6 Dashed red line is the ',size=1.8
			xyouts,40,0.4,/data,'!6 weighted mean PSF', size=1.8
			xyouts,40,0.36,/data,'!6 from Bookend iodines', size=1.8
			!p.color=0
			oplot,findgen(n_psf),vd_avpsf[i,*],col=222, thick=3, linesty=2
		  endif ;plt =2
	endfor ; i (all chunks0)
	   
; get the weighted mean wavelength solution (no smoothing) 
	for i=0,n_ord-1 do begin
		ordr=st_ord+i 
		ss=where(vdx.ordt eq ordr,nss)
 		str2='!6 Wavelength for Order '+strcompress(string(ordr),/rem)
 		
		for j=0,num_i2-1 do begin
		    st_ch_indx=j*nch_ord
		    end_ch_indx=(j+1)*nch_ord-1
	    ; print,'Chunk index across an order: ',ordr, st_ch_indx, end_ch_indx
	  	    w0=vdx[ss[st_ch_indx:end_ch_indx]].w0
		    wcof1=vdx[ss[st_ch_indx:end_ch_indx]].wcof[0]
		    wcof2=vdx[ss[st_ch_indx:end_ch_indx]].wcof[1]
		    wav0[i,*,j]=w0 + wcof1 ; + findgen(nch_ord)*wcof2
		    disp[i,*,j]=wcof2
		    dum = wt[ss[st_ch_indx:end_ch_indx],j]
	 	endfor  ; reading in wavelength solns for an order for each bookend
 		if plt gt 4 then begin
			xx2=findgen(nch_ord)
			for j=0,num_i2-1 do begin
				if j eq 0 then plot,xx2,wav0[i,*,0], xtitl=str2,/xsty,/ysty,ps=8,symsize=2.
				if j gt 0 then oplot,xx2,wav0[i,*,j],ps=8, col=j*45, symsize=2./j
			endfor
		endif	; plt = 4
	endfor ; each order

; Calculate the average vdiod chunk wavelength and dispersion
	chk_ct=0
	for m=0,n_ord-1 do begin
		for l=0,nch_ord-1 do begin
			vd_avwav[chk_ct]=median( wav0[m,l,*] )
			err_wav[chk_ct]=stddev(double(wav0[m,l,*]))
			vd_avdisp[chk_ct]= median( disp[m,l,*] )
			chk_ct+=1
		endfor 

		if keyword_set(smth) then begin
	; replace deviant vd_avwav and vd_avdisp values to modeled values
			xxarr=findgen(nch_ord)
			print,m*nch_ord, (m+1)*nch_ord-1
			ord_wav=vd_avwav[m*nch_ord:(m+1)*nch_ord-1]
			err_ord_wav=err_wav[m*nch_ord:(m+1)*nch_ord-1]>0.0001
		;remove dominant linear trend, then fit
			lin_coef=robust_poly_fit(xxarr,ord_wav,1,ylinfit,sig,/double)
			yresid=ord_wav-ylinfit
			nonlin_coef=robust_poly_fit(xxarr,yresid,6,ynlin,sig,/double)
			diff_wav=ord_wav-(ylinfit+ynlin)
			sbad_wav=where(abs(diff_wav) gt 0.005,nsbad_wav)
			if nsbad_wav gt 0 then begin
				ord_wav[sbad_wav] = ylinfit[sbad_wav]+ynlin[sbad_wav]	
;				ord_wav = ylinfit+ynlin
				print, 'Order: ',m,' number of adjusted wavelength points: ',nsbad_wav
				plot,xxarr,ord_wav,ps=8,/xsty,/ysty,$
						titl='!6 Wavelength zero points'
;			oplot,xxarr,ywfit,col=222
				oplot,xxarr[sbad_wav],ylinfit[sbad_wav]+ynlin[sbad_wav],ps=8, col=155
				wait,2
			;refit
				lin_coef2=robust_poly_fit(xxarr,ord_wav,1,ylinfit2,sig,/double)
				yresid2=ord_wav-ylinfit2
				nonlin_coef2=robust_poly_fit(xxarr,yresid2,6,ynlin2,sig,/double)
				diff_wav2=ord_wav-(ylinfit2+ynlin2)
				sbad_wav2=where(abs(diff_wav2) gt 0.001,nsbad_wav2)
			 	if nsbad_wav2 gt 0 then begin
					ord_wav[sbad_wav2] = ylinfit2[sbad_wav2]+ynlin2[sbad_wav2]	
					print,'second round bad wavelengths: ',nsbad_wav2
				endif
			endif ;fixing deviant wavelength points
			vd_avwav[m*nch_ord:(m+1)*nch_ord-1]=ord_wav

			if plt ge 3 then begin
				plot,xxarr,ord_wav,ps=8,/xsty,/ysty,$
							titl='!6 Wavelength zero points'
				oplot,xxarr,ylinfit+ynlin,col=222
				wait,1
			endif ; plt eq 2
		
		; now smooth the dispersion
		; first remove the dominent linear trend
			ord_disp=vd_avdisp[m*nch_ord:(m+1)*nch_ord-1]
			lin_disp_coef=robust_poly_fit(xxarr,ord_disp,1,ylindisp,sig,/double)
			ydisp_resid=ord_disp-ylindisp
			; then fit the residual nonlinear curve
			nl_disp_coef=robust_poly_fit(xxarr,ydisp_resid,4,ydisp_nlin,sig,/double)
			diff_disp=ord_disp-(ylindisp+ydisp_nlin)
			sbad_disp=where(abs(diff_disp) gt 0.0005,nsbad_disp) 
			if nsbad_disp gt 0 then begin
				ord_disp[sbad_disp] = ylindisp[sbad_disp]+ydisp_nlin[sbad_disp]
;				ord_disp = ylindisp + ydisp_nlin
				print, 'Order: ',m,' number of adjusted dispersion points: ',nsbad_disp
				plot,xxarr,ord_disp,ps=8,col=0,/xsty,/ysty,$
					titl='!6 Dispersion'
;				oplot,xxarr,ydisp,col=222
				oplot,xxarr[sbad_disp],ylindisp[sbad_disp]+ydisp_nlin[sbad_disp],ps=8,col=155
				wait,2
				;refit
				lin_disp_coef2=robust_poly_fit(xxarr,ord_disp,1,ylindisp2,sig,/double)
				ydisp_resid2=ord_disp-ylindisp2
				; then fit the residual nonlinear curve
				nl_disp_coef2=robust_poly_fit(xxarr,ydisp_resid2,4,ydisp_nlin2,sig,/double)
				diff_disp2=ord_disp-(ylindisp2+ydisp_nlin2)
				sbad_disp2=where(abs(diff_disp2) gt 0.0005,nsbad_disp2) 
			 	if nsbad_disp2 gt 0 then begin
					ord_disp[sbad_disp2] = ylindisp2[sbad_disp2]+ydisp_nlin2[sbad_disp2]	
					print,'second round bad dispersion: ',nsbad_disp2
				endif
			endif ;fixing deviant dispersion points		
			vd_avdisp[m*nch_ord:(m+1)*nch_ord-1]=ord_disp
			if plt ge 3 then begin
				plot,xxarr,ord_disp,ps=8,col=0,/xsty,/ysty,$
					titl='!6 Corrected dispersion'
				oplot,xxarr,ylindisp+ydisp_nlin,col=222
				wait,1
			endif ; plt eq 3
		endif ; keyword_set(smooth)
	endfor  ;m (looping through orders)

endif  ;not vdavg

if keyword_set(vdavg) then begin
	yrmo=strmid(date,0,4)
	restore,'/tous/mir7/files_df/vd'+vdtag+'AVG'+'_'+yrmo  ;vdav
	n_psf=n_elements(vdav[0].psf)      	; number of array elements in the psf (typically 121)
	vd_avpsf=dblarr(n_chunks, n_psf)  	; weighted ave psf for each chunk
	vd_avwav = dblarr(n_chunks)       	; weighted mean wavelength soln for each chunk 
	vd_avdisp = dblarr(n_chunks)      	; weighted mean dispersion soln for each chunk 
	vdx=vdav 
	for i=0,n_chunks-1 do begin
		vd_avpsf[i,*]=vdav[i].psf
		vd_avwav[i]=vdav[i].w0 + vdav[i].wcof[0]
		vd_avdisp[i]=vdav[i].wcof[1]
	endfor  ;nchunks
	num_i2=10
endif  ;vdavg chunks
	
dsst_nm='dsst'+iss_nm+iss_tag+'_a'+obsnm+'.dat'
print, ' Analyzing '+dsst_nm 

restore,dopcode_path+pfilt_file  ;filt 

; need some pixel padding for the deconvolution
if ~keyword_set(dsst_pad) then dsst_pad=12   ; extra pixels
dstlen = osamp*(npix + 2*dsst_pad)
tmp_pad=30  ; temporary padding for continuum normalization
if tmp_pad lt dsst_pad then stop,'dsst_pad is too large' 

; dsst structure, fill with ordt and pixt values
dsst={ordt:0,pixt:0,w0:double(0.),w1:double(0.),weight:0.d,dst:dblarr(dstlen)}
dsst=replicate(dsst,n_chunks) 
dsst.ordt = vdx[0:n_chunks-1].ordt
dsst.pixt = vdx[0:n_chunks-1].pixt

; fill dsst with wavelength soln, weights and deconvolved spectrum (dst)
for i=0,n_chunks-1 do begin  ; fill the dsst structure
	;find w0 for the padded dsst chunk
	lambda = vd_avwav[i]       ; wav for pix0 (not padded)
	dispersion = vd_avdisp[i]  ; disp for pix0
	ord=dsst[i].ordt
	pix0=dsst[i].pixt

	; now get coadded template observation
	restore,dfd+coadded_spec  ;star array
	; calculate weight from line slopes
	; (G Marcy analytical weights) single value for dsst chunk
	sp=star[pix0:pix0+npix-1, ord]
	eps=sqrt(sp)
	didp=sp[1:npix-1] - sp[0:npix-2]    ; slope dI/d(pix)
	didv=didp*lambda/(c*dispersion) ; slope in Intensity per m/s
	dsst[i].weight=double(total( (didv/eps)^2 ))

    ; fill the dsst.dst with the spectrum 
    ; must send deconv_j a normalized spectrum!
    ; first do a rough continuum normalization with padded chunks
    sz1=size(star)  
    maxpix=sz1[1]-1       ;3200 pixels per order
    minpix=0.
    lo=pix0-tmp_pad > minpix     ; 50
    hi=pix0+(npix-1)+tmp_pad < maxpix  ;189
    pix_filt=filt[lo:hi,ord]
 	nf=n_elements(pix_filt)
 	    
	; zero weight for the edges of the chunk
 	pix_filt[0:2]=0
 	pix_filt[nf-3:nf-1]=0

 	strr=reform(star[*,ord])
	flstrr=strr
;	llo=0 > (pix0-tmp_pad)     ; 200 elements, [20, 219]
;	hhi=(llo+((npad/2)*tmp_pad)+(npix-1) + ((npad/2)*tmp_pad)) < maxpix
	tmp_spec = strr[lo:hi]           ; 140 elements 

	; continuum normalization
    contf,tmp_spec,c_tmp,sbin=20,nord=1, frac=0.1
    tmp_spec=tmp_spec/c_tmp
    c2=max(tmp_spec)
    cont_scale = mean(c2)/mean(c_tmp)
    spec=tmp_spec/c2               ; 140 elements [50, 189]
;    spec = spec[(tmp_pad-dsst_pad):(n_elements(spec) - 2.*(tmp_pad - dsst_pad) -1)]
	nspec=n_elements(spec) 
	     
	wdum=lambda+dispersion*(indgen(nspec)-tmp_pad)  ;wav for spec and tmp_spec
	wv0=lambda-dispersion*dsst_pad                  ;wav for dsst.dst
    psf=reform(vd_avpsf[i,*])

 	if keyword_set(jansson) then begin
 	; spec has tmp_pad on each end (more than dsst_pad), which will return 
 	; 80 + 2.*30 = 140 pixels or 560 osamp pixels
 	; after deconvolution, trim this down to 80 + 2.*dsst_pad = 120 pixels or 480 osamp pixels
 	
	;   This routine deconvolves a spectrum (spec)
	;     using the Jansson constrained non-linear
	;     method (Deconvolution With Applicatons
	;     in Spectroscopy, by P.A. Jansson,
	;     Academic Press, 1984, QC451.6, .D45),
	;     as outlined by Ronald L. Gilliland,
	;     Space Telescope, Baltimore.
	;
	;   Written  by Paul Butler, July 22, 1991
	;   Modified by Paul Butler, July 20, 1993 to allow for "parabolic smoothing"
	;   Modified by Paul Butler, Sept 24, 1995 to allow for "direct oversampling"
	;  Oct 2001, DAF - idl's internal fspline is now being called instead of 
	;  Jul 2009, DAF - revised to work with new Doppler code and CTIO
	;                  spectral format
	; Jan 2013, DAF - revised to use the average PSF from yrmo Bstars
	
   ; smooth the spectrum to reduce ringing in the continuum
   	sm_spec=smooth(spec, smooth_val)
   	;print,stddev((spec - sm_spec)/spec)
   	if stddev((spec - sm_spec)/spec) gt 0.5 then print,'smoothing failed' 
   	spec=sm_spec

 		parsm=6 ;parabolic smoothing for Jansson deconvolution
		dsst_nm='dsst'+iss_nm+iss_tag+'jansson_a'+obsnm+'.dat'

		nspec = n_elements(spec)                ;length of input spectrum
		xx  = findgen(nspec)                    ;useful array, x-array for spec
		sxx = xx  ;[1:nspec-2]                     ;x-array for spec without endpoints
		xxx = findgen(nspec*osamp)/float(osamp)  ;x-array for oversampled result
		sp  = reform(float(fspline(xx,spec,xxx)))  ;first guess for oversampled deconvlved spec
		if n_elements(pix_filt) ne n_elements(spec) then pix_filt=spec*0+1
		if n_elements(j_bval) ne 1 then b=1. else b=j_bval
		a=0.    &     aplusb=a+b   &   bminusa=b-a  
		rk=1.                                   ;Corrections application, see Jansson, Gilliand, etc ....
		if n_elements(parsm) ne 1 then parsm=1  ;default, no parabolic smoothing of difference array
		if parsm eq -1 then print,'Deconv_j: Fourier smoothing of difference array'
		if parsm gt 1 then begin                ;differnce array is parabolicly smoothed
	    parsm=fix(abs(parsm))
	    nosxx=2*parsm+1                      ;Number of OverSampled piXXXels
   		osxx=findgen(nosxx) - parsm          ;OverSampled X array 
   		mxx=fltarr(3,nosxx)                  ;Parabolic X array Matrix (linear algebra)
   		for m=0,2 do mxx[m,*]=osxx^m         ;Setting up Parabolic X array
		   t_xarr=transpose(mxx)                ;Linear algebra manipulations ...
		   sqmat=mxx # t_xarr                   ;Linear algebra manipulations ...
		   invsqmat=invert(sqmat)               ;Linear algebra manipulations ...
		   mxx=t_xarr # invsqmat                ;Parabolic smoothing matrix
		endif   ;parsm

		nln=n_elements(sp)                      ;number of oversampled pixels
		ln=long(nln)-1
		dsp=sp  &     fake_sp=sp

		;*** Convergence Criteria ***
		chi_change=.001  &   ctmax=100   

		;  strip instrumental profile down to "non-zero" region
		;  and do preliminary convolution setup
		inpr=psf
		ind=where(abs(inpr) gt 0.001*max(inpr),n_ind)
		lo_ip=max([ind[0]-1,0])  &  hi_ip=min([ind(n_ind-1)+1,n_elements(inpr)-1])
		inpr=inpr[lo_ip:hi_ip]
		inpr=reverse(inpr)/total(inpr)

		fltind=where(pix_filt gt 0)           ;filter index
		r_conv,dsp,inpr,fake_sp               ;convolving first guess deconvolved oversampled spectrum
		rebin,xxx,fake_sp,sxx,fake_spec       ;binning to standard sampling
        ndiff=fltarr(nspec)
		ndiff[fltind]=spec[fltind]-fake_spec[fltind]  ;difference array, no oversampling
		diff = fspline(xx,ndiff,xxx)          ;oversampling the difference array
	
		;  set up iteration loop
		chi=1.
		icount=-1  ; iteration counter
		repeat begin
                   if parsm gt 1 then begin ;parabolic smoothing differnces
                      dumdiff=diff
                      for qq=parsm,(ln-parsm) do $ ;Parabolic smoothing!   
                         dumdiff(qq)=first_el(reform(diff(osxx+qq)) # mxx ) 
                      diff=dumdiff
                   endif
   
                 ;print,icount
                  icount=icount+1    &   old_chi=chi

		 ;Generate the next iteration of the deconvolved oversampled spectrum
                  rcor=rk*(1-abs(dsp-aplusb/2.)*2./bminusa) ;corrections
                  dspec=dsp+rcor*diff  &  dsp=dspec         ;new deconvolved oversampled spectrum
                  r_conv,dsp,inpr,fake_sp                   ;convolution of oversampled deconvolved spectrum
                  rebin,xxx,fake_sp,sxx,fake_spec           ;binning to standard sampling
                  ndiff=spec-fake_spec                     ;differences in standard sampling space
                  diff = fspline(xx,ndiff,xxx)             ;differences in oversampled space
                  chi=stddev(ndiff(fltind))                ;standard deviation
                  ;  Stop iterating when EITHER fractional change in chi lt chichange
                  ;  or icount reaches ctmax

		endrep until ((old_chi/chi-1.) lt chi_change or icount ge ctmax)
		dspec=dsp                             ;final deconvolved oversampled result
		if plt ge 1 then begin
			plot,sp,linestyle=0,thick=2
;			oplot,rcor,col=222
			oplot,dspec,col=80,thick=1
		endif 

		wstr=findgen(n_elements(dspec))*(dispersion/float(osamp)) + min(wdum)
		dd=double(abs(wstr - wv0))
		l=where(dd eq min(dd),nl)
		px0 = l
		dspec=dspec[px0:px0+dstlen-1]
		wstr=wstr[px0:px0+dstlen-1]
		dsst[i].dst = dspec
		dsst[i].w0 = wstr[0]
		dsst[i].w1 = dispersion/float(osamp) 
		wfine = dsst[i].w0 + dsst[i].w1*indgen(dstlen)

endif ;jansson

	if ~keyword_set(fun) and ~keyword_set(jansson) then begin  ;use mpfit (not mpfitfun) for better results
        ; a large number of nodes introduces too much structure
        ; assume 2-pix sampling for CHIRON's observed spectrum
		if ~keyword_set(ngrid) then ngrid=20  ;(number of nodes for dop_deconv_df) 

        xfine = makearr(nspec*osamp, 0, nspec-1./osamp)  ; osamp pixel steps
		osamp_spec = dspline(indgen(nspec), spec, xfine)
		nfine = n_elements(xfine) 
		
		if 1-keyword_set(sigma) then sigma = 3.0 ;;; Characteristic width (pixels)
		sig = 1 ;;; pixel unc. 
        
		; no nodes at the edges of the chunk
		if 1-keyword_set(lim) then lim = [5, float(nfine)/float(osamp)-5]
		node = makearr(ngrid, lim[0], lim[1])
		; to fill dsst.dst:
		; need to generate the function dspec (deconvolved spectrum) such that 
		; dspec convolved with the PSF => observed template spec

		par = dblarr(2*ngrid+1)  ; each node has amplitude and offset as free parameter
		par[0] = sigma            ; line resolution
		npar = n_elements(par)    
		trim = 10

		  fa = {xfine:xfine,       $
				spec:osamp_spec,   $
      			psf:psf,           $
      			sig:sig,           $
				node:node,         $
				trim: trim,        $
				mov: mov,          $
				ngrid: ngrid}
     			
		newpar = mpfit('dop_deconv_df', par, functargs=fa, parinfo=parinfo,$
					/quiet, maxiter=300)
		dum=dop_deconv_df(newpar, xfine=fa.xfine, spec=fa.spec, psf=fa.psf, sig=fa.sig,$
				node=fa.node, ngrid=fa.ngrid, trim=fa.trim, newspec=dspec, resid=resid)
				;dspec is the deconvolved spectrum

		sp=star[pix0:pix0+npix-1, ord]		
		if plt eq 1 then begin
			plot,sp,linestyle=0,thick=2
;			oplot,rcor,col=222
			oplot,dspec,col=80,thick=1
		endif 

		wstr=findgen(n_elements(dspec))*(dispersion/float(osamp)) + min(wdum)
		dd=double(abs(wstr - wv0))
		l=where(dd eq min(dd),nl)
		px0 = l
		dspec=dspec[px0:px0+dstlen-1]
		wstr=wstr[px0:px0+dstlen-1]
		dsst[i].dst = dspec
		dsst[i].w0 = wstr[0]
		dsst[i].w1 = dispersion/float(osamp) 
		wfine = dsst[i].w0 + dsst[i].w1*indgen(dstlen)		
	endif

;	if keyword_set(fun) then begin  ;not using this section now.
;		wt=1./(osamp_spec*cont_scale)
;		wt=wt/total(wt) 
;		newpar = mpfitfun('dop_deconv_dfun', xfine, osamp_spec, parinfo=parinfo, $
;				functargs=fa, maxiter=100, $
;				yfit=dspec,weight=wt, status=status )
;		dum = dop_deconv_dfun(xfine, newpar, functargs=fa)
;		; having trouble getting dspec returned from mpfitfun! help? 
;		; no "fun" for now... 
;	endif

;	if keyword_set(movie) then wset,0

		if plt ge 1 then begin
        	xt='Wavelength [A]'
        	yt='Intensity'
        	ttl='Order '+strtrim(ord,2)+', Pixel '+strtrim(pix0,2)
        	plot,wdum, spec, /xsty,/ynoz, xtitle=xt,ytitle=yt,title=ttl
        	oplot,wstr, dspec, col=220
        	dd=convol(dspec, psf, /edge_truncate,/normalize)
        	conv_spec=dspline(wstr, dd, wdum)
			xgd=where(wdum ge wstr[0] and wdum le wstr[n_elements(wstr)-1] )
        	oplot,wdum[xgd],conv_spec[xgd],col=155,ps=-8,symsize=0.5
        	wait,0.1
    	endif
    	
    if i eq 0 and ~keyword_set(noprint) then begin
        print,' '
        print,'                  *  CREATING DSST  *     '
        print,' '
        print,'|--------------------------------------------------------|'
        print,'| Order Pixel  Lambda   dlam/dx  Scat    SLOPE  # PSFs   |'
        print,'|               (A)     (A/pix)           SUM   avg''d    |'
        print,'|--------------------------------------------------------|'
    endif
    if keyword_set(noprint) then begin
        counter, n+1, n_lines, 'Chunk # ',/timeleft,starttime=stt
    endif else begin
        fmt = '(I5,I6,F10.3,F9.5,F8.3,F8.2,I6)'
        print,format=fmt,ord,pix0,dsst[i].w0,dsst[i].w1,0,dsst[i].weight,num_i2
    endelse

endfor  ;chunk loop

print,dsst_nm
save,dsst, f=dsst_nm

if plt eq 1 then begin  ;is the dsst continuous?
cc=40
	wd1=double(dsst[cc].w0 + indgen(dstlen)*dsst[cc].w1)
	wd2=double(dsst[cc+1].w0 + indgen(dstlen)*dsst[cc+1].w1)
	xr1=min(wd1)  & xr2=max(wd2)
	plot, wd1, dsst[cc].dst, /ysty, xra=[xr1,xr2],/xsty
	oplot, wd2, dsst[cc+1].dst, col=222
endif

ans='n'
read,'save to /tous/mir7/files_df/? (y/n): ',ans
if ans eq 'y' then save, dsst, f='/tous/mir7/files_df/'+dsst_nm
stop
end

