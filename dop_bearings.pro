;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PRO DOP_BEARINGS
;
; PROCEDURES CALLED: MY_INIT, (e.g., CTIO_INIT), DOP_INIT, DOP_LOC, PSF
;                    Also need wavelength soln, 1 wavelength assigned
;                    to each pixel in the 2d array passed in "wavfile"
;
; PURPOSE: 
;   To help you get the wavelength solution for starting a new program
;
; PROCEDURE:
;   1) gets the set up (pathnames and input variables from your version of 
;      MY_INIT, e.g., CTIO_INIT) 
;   2) loads dopenv, and IDL structure with program specific parameters
;   3) compares the observed I2 with a (convolved) FTS wavelength segment
;
; CALLING SEQUENCE: 
;   dop_bearings, obsnm='rqa06.7452', wavfile='ctio_rqa06.7400.dat' 
;
; INPUTS: 
;   OBSNM:  (obs name) stellar spectrum to be analyzed
;   WAVFILE: wavelength soln, 1 wavelength for each pixel in the 2-d
;   array 
;   TAG: to distinguish one trial run from another, a tag is added to
;   the output files - the tag is arbitrary, and a string. 
;   INP_VD: untested - the idea was to make a second pass, but
;                      doesnt' seem to be necessary.  ignore for now. 
;
; OUTPUTS
;   VD file: eg "vdciod_rg05.1351" this structure contains the basic
;   wavelength soln 
;   CHUNK file: a structure with configured with just the barebones
;   information in the VD file, but the compatible match for the
;   doppler code. 
;
; Written by Debra Fischer, SFSU, Nov 2009
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO DOP_BEARINGS, obsnm=obsnm, wavfile=wavfile, tag=tag, observatory=observatory,$
     inp_vd=inp_vd,delta=delta, mode=mode, rev_thid=rev_thid, wid=wid, date=date

; fischer oct 2009
; set up keck vd structure for the first time
; this code might bypass the need for an ipcf file

; INPUT
;   obsnm: Bstar iodine observation
;   wavfile, wavelength soln, e.g. ctio_rqa06.7400.dat
;   OPTIONAL: inp_vd (refine existing solution) - NOT TESTED
;
;   dop_bearings,obsnm='ri67.220',wavfile='/home/fischer/thid/lick_ri67.214.dat', tag='e',observ='lick'
;

if mode eq 'slit' then wid=1.95
if mode eq 'narrow_slit' then wid=1.5

if ~keyword_set(observatory) then begin
   observatory=''
   read,'which observatory (lick, keck, het, ctio, ctio4k): ',observatory
endif

     if observatory eq 'ctio4k' then dopenv=ctio4k_init(obsnm, 'iod', 0.0, mode=mode, iss_obnm=obsnm,tag=tag,date=date) 
     if observatory eq 'ctio' then dopenv=ctio_init(obsnm, 'iod', 0.0, iss_obnm=obsnm,tag=tag) 
     if observatory eq 'keck' then dopenv=keck_init(obsnm,  'iod', 0.0, iss_obnm=obsnm,tag=tag) 
     if observatory eq 'lick' then dopenv=lick_init(obsnm,  'iod', 0.0, iss_obnm=obsnm,tag=tag) 
;     if observatory eq 'het' then dopenv=het_init(obsnm, 'iod', 0.0, iss_obnm=obsnm,tag=tag) 
     nchunks=dopenv.n_chunks
     if ~keyword_set(tag) then begin
        tag=''
        read,'which tag (e.g. b): ',tag
     endif

     npsf=n_elements(dopenv.psfpix)
     vdst={obnm:dopenv.obsnm, ordt:0, ordob:0, pixt:0, pixob:0, w0:0, wcof:dblarr(dopenv.n_wcof), cts:0L, $
           z:0.0d, fit:0.0, npix:160, vel:0.0, weight:0.0, $
           psf:fltarr(121), par:fltarr(npsf)}
     vd=replicate(vdst, nchunks)
    
     nord=dopenv.n_ord                    ; number of useable iodine orders 
     npix=dopenv.n_pix                    ; number of pixels per chunk

     if ~keyword_set(het) then bst=mrdfits(dopenv.obs_dir+obsnm+'.fits',dopenv.file_ext) 
     if keyword_set(het) then bst=mrdfits(dopenv.obs_dir+obsnm+'.fits.gz',dopenv.file_ext)
     if n_elements(size(bst)) gt 5 then bst=reform(bst[1,*,*])  ;wavelength array was attached
      ;bst=reverse(bst,2)
                                ;stop  ; is the keck fits spectrum
                                ;reversed?  yes, if fits files are
                                ;from berkeley
     restore,wavfile          ; restore wav array "w" 
     if keyword_set(rev_thid) then w=reverse(w,2)
     if ~keyword_set(delta) then read,'enter delta (e.g., 1.4): ',delta
;     if observatory eq 'ctio' then w=w*(1.+delta/w[0,5])
     w=w*(1.+delta/w[0,5])   ;applied to the entire wavelength soln as a shift
     jdisp = (w[dopenv.st_pix+dopenv.n_pix,dopenv.st_ord] - w[dopenv.st_pix,dopenv.st_ord])/dopenv.n_pix

   ; THE FOLLOWING BLOCK IS NOT YET TESTED
   ; if an input vd exists (already run this code) then 
   ; you can fine tune the wavelength solution with a second pass
   ;  if keyword_set(inp_vd) then begin
   ;     restore,inp_vd
   ;     iuniq=uniq(vd.ordt)
   ;     ord=vd[iuniq].ordt
   ;     if n_elements(ord) ne nord then stop,'why not the same??'
   ;     for ii = 0,n_elements(ord)-1 do begin
   ;        if ii eq 0 then begin
   ;           xch=where(vd.ordt eq ii,nxch) ;12 of these
   ;           for jj=0,nxch-1 do begin
   ;              wav=fltarr(npix)
   ;              stw=950+jj*npix
   ;              w[stw:stw+npix,ii]=(vd[jj].w0 + vd[jj].wcof[0]) + vd[jj].wcof[1]*findgen(npix)
   ;           endfor ;loop for first order              
   ;        endif ;ii = 1
   ;        if ii gt 0 then begin
   ;           xch=where(vd.ordt eq ii,nxch) ;35 of these
   ;           for jj=0,nxch-1 do begin
   ;              wav=fltarr(npix)
   ;              stw=950+jj*npix
   ;              w[stw:stw+npix,ii]=(vd[jj].w0 + vd[jj].wcof[0]) + vd[jj].wcof[1]*findgen(npix)
   ;           endfor ;loop for next orders              
   ;        endif ;ii gt 1
   ;     endfor  ; done looping through the ordt's
   ;stop  ; check the wavelength array
   ; endif


    ; FIRST PASS - MATCH THE BSTAR IOD SPEC TO CONVOLVED FTS SPEC
    ; THIS LOOP CHECKS SOLN AND "PATCH" FOR EACH CHUNK
    ; FILLS WAVE[*,J] AND CHUNK[*,J]
    kcount=0
    !p.font=0
    !p.charsize=2
    delta2 = 0.0      ; shift for wavelength offset, determined for each order
    count_ord=0
    restore,'torrent_filt.dat'

  for i=dopenv.st_ord, dopenv.st_ord+nord - 1 do begin   ;loop through each of the orders
       wtmp=reform(w[*,i])      ; one order at a time
       btmp=reform(bst[*,i])    ; one order at a time
       ;stop
       st_pix=dopenv.st_pix
       nch_ord=dopenv.nch_ord
       if observatory eq 'keck' and i eq 0 then begin
          st_pix=1410
          nch_ord=32
       endif
       if observatory eq 'keck' and i gt 0 then begin
          st_pix=50
          nch_ord=49
       endif
       if observatory eq 'keck' and i eq 14 then begin
          st_pix=50
          nch_ord=33
       endif

       if observatory eq 'het' and i eq 2 then begin
          st_pix=990
          nch_ord=29
       endif
       if observatory eq 'het' and i eq 2 then begin
          st_pix=990
          nch_ord=30
       endif

       if observatory eq 'het' and i eq 19 then begin
          st_pix=1980
          nch_ord=21
       endif

loadct,39
       for j=0,nch_ord-1 do begin   ; loop through chunks in the order
          istart=j*npix + st_pix
          iend=istart+(npix-1)
          wtmp_ch=wtmp[istart:iend] ;thid wavelength
          btmp_ch=btmp[istart:iend]
          cts=median(btmp_ch)
          pfilt=filt[istart:iend, i]
          noise_flat=0.008
          wt=(1./btmp_ch) / (1. + btmp_ch*noise_flat^2)
          wt = wt*pfilt*1.0d    ; zero out the bad pixels 
          sig=1./sqrt(wt)
          contf,btmp_ch,cb,sbin=20,nord=1 ;just 80 pixels so nord=1 is OK
          btmp_ch=btmp_ch/cb

          if ~keyword_set(wid) then wid = 2.3            ; guess the FWHM for the convolution 

        ; get the FTS spectrum for that same wavelength interval
          if observatory eq 'ctio4k' then read_iodine, 'ctiov1', 'pnnl', 40, min(wtmp_ch), max(wtmp_ch), wiod, siod, nadd=50
          if dopenv.fts_atlas eq 'lickq1_pnnl_red10.sav' then begin
             read_iodine,'lickq1','pnnl',50.0, min(wtmp_ch), max(wtmp_ch),wiod,siod,nadd=100 
          endif 
;             rdfts,wiod,siod,min(wtmp_ch), max(wtmp_ch),dfd=dopenv.fts_dir,dfn=dopenv.fts_atlas
            if dopenv.fts_atlas eq 'iodine_p4_pnnl_lo_apod.dat' or $  ; keck
            	dopenv.fts_atlas eq 'iodine_p4_pnnl_lolo.dat' then $
		  		read_iod_p4,dopenv.fts_atlas,wiod, siod, min(wtmp_ch), max(wtmp_ch),pad=1.
             contf,siod,c,sbin=20,nord=1
             siod=siod/c
          disp=(wtmp_ch[npix-1] - wtmp_ch[0])/npix
          shift=0.0             ; think about adding this with cursor command input
          offset=0.99
          par=[wtmp_ch[0], disp, shift,1.0]
          npar=n_elements(par)
           
          if i eq dopenv.st_ord and j eq 0 then help = 1 else help=0
          ans='n'  &   delta=0.
          if help eq 1 then begin
             while ans eq 'n' do begin 
                                ; quick look check - maybe make this
                                ;                    just for the
                                ;                    first chunk in
                                ;                    each order? 
                print,'HELP ME GET THE WAVELENGTH SOLUTION' 
                xarr=findgen(npix) ; number of pixels for rebinning the siod
                ip=psf(wid)
                testspec=convol(siod,ip,/edge_truncate,/normalize)
                plot,wiod,testspec-0.2,/xsty,xtitl='!6FTS wavelength [A]',ytit='!6 Convolved FTS Spec',yra=[0.5,1.3]
                                ;           oplot,wiod,siod,col=1,linesty=2
                oplot,wtmp_ch,btmp_ch,col=222,thick=2
                !p.color=1
                xyouts,max(wtmp_ch)-1.0,1.2,/data,'!6 Black is FTS',size=2.
                !p.color=220
                xyouts,max(wtmp_ch)-1.0,1.1,/data,'!6 Red is obs',size=2.
                !p.color=1

                print,'Click on any line in the Observation (red)'
                cursor,wtemp,y
                wait,1
                print,'Click on the corresponding line in the NSO (black)'
                cursor,wnso,y
                wait,1                
                delta2=wnso - wtemp
                wtmp_ch=wtmp_ch*(1. + delta2/wtmp_ch[0])
                plot,wiod,testspec-0.1,/xsty,xtitl='!6FTS wavelength [A]',ytit='!6 Convolved FTS Spec',yra=[0.5,1.3]
                                ;           oplot,wiod,siod,col=1,linesty=2
                oplot,wtmp_ch,btmp_ch,col=222,thick=2
                read,'Is this the right shift? (y/n) ',ans
            end

             par[0] = par[0]*(1.+delta2/wtmp_ch[0])
             par[2] = 0.0 

             print,'NOW, HELP ME GET A GUESS FOR THE DISPERSION BY CLICKING ON TWO LINES'
             print,'Click on any blueward line in the Observation (red)'
             cursor,wtemp_bl,y
             wait,1
             diff1=abs(wtmp_ch-wtemp_bl)
             ind1=where(diff1 eq min(diff1))
             print,'Click on the corresponding blueward line in the NSO (black)'
             cursor,wnso_bl,y
             wait,1
             print,'Click on any redward line in the Observation (red)'
             cursor,wtemp_red,y
             wait,1
             diff2=abs(wtmp_ch-wtemp_red)
             ind2=where(diff2 eq min(diff2))
             npix_disp=ind2-ind1
             print,'Click on the corresponding redward line in the NSO (black)'
             cursor,wnso_red,y
             wait,1
             disp = (wnso_red - wnso_bl)/npix_disp
             par[1] = disp[0]
             wtmp_ch_new=par[0] + par[1]*findgen(npix)
             plot,wiod,testspec,col=1,charsize=1.8
             oplot,wtmp_ch_new,btmp_ch,col=222
          endif                 ; first chunk in the order gives delta w0 for next chunks

          functargs={btmp_ch:btmp_ch, wiod:wiod, siod:siod, wid:wid} ;structure of extra's for mpfit
          parinfo={value:0.0d,  $
                   parname:'?', $
                   step:0.01, $
                   fixed:0,   $
                   limited:[0,0], $
                   limits:fltarr(2) }
          parinfo=replicate(parinfo, npar)
          parinfo[0].value=par[0]   &  parinfo[0].parname='w0'  & parinfo[0].step=0.001
          parinfo[1].value=par[1]   &  parinfo[1].parname='disp'  & parinfo[1].step=0.0001
          parinfo[2].value=par[2]   &  parinfo[2].parname='wavelength shift' & parinfo[2].fixed=1
          parinfo[3].value=par[3]   &  parinfo[3].parname='offset'  & parinfo[3].step=0.1

          x=findgen(npix) 
          y=btmp_ch

          newpar=mpfitfun('dop_loc', x, y, parinfo=parinfo, $
                          functargs=functargs, weight=wt, yfit=syn_fit,/verbose) 
          wav=newpar[0] + newpar[1]*(x+newpar[2])
          plot,wav,syn_fit,/xsty
          oplot,wav,btmp_ch,col=222 
          !p.color=1
          xyouts,max(wav)-1.,0.45,/data,'!6 Convolved FTS Spectrum',size=1.8
          !p.color=222
          xyouts,max(wav)-1.,0.38,/data,'!6 Observed Spectrum',size=1.8
          !p.color=1
          red_chi_sq = total(cts^2*wt*(syn_fit - y)^2)/(npix-npar)
          print,'Reduced chisq: ',sqrt(red_chi_sq)
          
        ; NOW BUILD VST TAGS FOR A GIVEN ORDER 
          vd[kcount].obnm=obsnm
          vd[kcount].ordt=i
          vd[kcount].ordob=i
          vd[kcount].pixt=istart 
          vd[kcount].pixob=istart
          vd[kcount].w0=fix(newpar[0])
          vd[kcount].wcof[0]=newpar[0]-vd[kcount].w0
          vd[kcount].wcof[1]=newpar[1]
          vd[kcount].cts=cts
          vd[kcount].fit=sqrt(red_chi_sq)
          vd[kcount].par=dopenv.psf_ampl
          print,vd[kcount].w0, vd[kcount].wcof[0], vd[kcount].wcof[1], vd[kcount].pixt
          kcount=kcount+1  ; increment the chunk counter
       endfor                   ; chunk loop, nch_ord
        xvd=findgen(nch_ord)         
        stch= kcount - nch_ord
        endch= kcount -1
        count_ord=count_ord+1    ;increment the orders
        ;dispersion should be continuous for these first guesses
        disp=vd[stch:endch].wcof[1]
        coef=poly_fit(xvd,disp,1)
        plot,xvd,disp,ps=4,col=1,symsize=1.5
        oplot,xvd,poly(xvd,coef),ps=8,col=222
        vd[stch:endch].wcof[1]=poly(xvd,coef)

        ;wavelengths should be continuous for these first guesses
        w0=vd[stch:endch].wcof[0]+vd[stch:endch].w0
        pix=vd[stch:endch].pixt
        lcoef=poly_fit(pix,w0,1,/double,yfit=lin_lpix)
        nlwav = w0 - lin_lpix
        poly_ord=6
        resid_lcoef=poly_fit(w0,nlwav,poly_ord,/double,yfit=lwavfit)
        vd[stch:endch].w0 = fix(double(lwavfit+lin_lpix))
        vd[stch:endch].wcof[0] = double(lwavfit+lin_lpix)-vd[stch:endch].w0
     endfor                     ; order loop, i 

  vdiodnm='vd'+tag+dopenv.obj_nm+'_'+obsnm
  print,vdiodnm
  save,vd,f=dopenv.files_dir+vdiodnm 
  conv_vd2cb,dopenv, tag=tag, vdiodnm=vdiodnm, obsnm=obsnm

end

