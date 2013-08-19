;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_MARQ
;
; CALLED FROM: DOP_MAIN
;
; CALLS TO: DOP_FIT (FIRST PASS) 
; 
; PURPOSE:  
;   Drives Levenberg-Marquardt search for best fit parameters:
;   PSF: Point Spread Function 
;   ISS: Intrinsic Stellar Spectrum
;   FTS I2: high resolution FTS spectrum of I2 cell
;
;   ISS * FTS_I2 convolved with PSF to model program observation
;    
; PROCEDURE:
;   1) weight the pixels 
;   2) set up the parinfo structure for free parameters
;   3) pre-pass: wavelength soln, doppler shift, 
;      width of central gaussian
;   4) first pass: everything free: wavelengths, doppler, PSF 
;   5) second pass: average the PSF
; 
; INPUTS: 
;      sobs: observed spectrum
;      siod: FTS spectrum;        wiod: FTS wavelengths
;      sis: intrinsic spectrum;  wis: wavelengths of intrinsic
;      psfpix, psfsig (position and widths of Gaussian's for PSF
;      chunk.free_pars
; 
; OUTPUTS:
;   MODEL: chunk.free_pars (updated) 
;   
; OUTSTANDING ISSUES: 
;    freeze dispersion on a later pass
;
; Written by: Debra Fischer, SFSU, Dec 2007
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_MARQ, chunk, dopenv, pass=pass, ch_indx, chunk_arr=chunk_arr,$
                   demo=demo, verbose=verbose, lm=lm, pfilt=pfilt, avg=avg

  zp = 0  ; zeroth Doppler pass par's
  fp = 1  ; first Doppler pass par's
  lp = 2  ; last Doppler pass par's

  ;PERTURBATIONS:
  ;Set this to 1 to slightly perturb the initial guesses. Otherwise, set to 0
  perturb = 1

    ; PIXEL WEIGHTS FOR CALCULATING CHISQ OF THE MODEL	
          ;noise_flat = 0.008    ;typical s/n in flatfield
          noise_flat = 0.005    ;typical s/n in flatfield (DF jun30 2012) 
          wt = (1./(*chunk.sobs)) / (1. + (*chunk.sobs)*noise_flat^2)
          xneg=where(wt lt 0.,nxneg) & if nxneg gt 0 then wt[xneg]=0.d
          xhi=where(wt gt 5.*median(wt),nxhi) & if nxhi gt 0 then wt[xhi]=0.d 
          wt = wt*pfilt*1.0d   ; zero out the bad pixels 
		  sig = sqrt(*chunk.sobs) ;only used in plotting

    ; BEGIN L-M FITTING
    ; INFORMATION AND LIMITS FOR FREE PARAMETERS: PARINFO STRUCTURE

    if pass eq 1 then begin   ;**** BEGIN PASS 1 ****

	if dopenv.n_wcof eq 3 then ij=1 else ij=0
        
       ; ******* QUICK REFINING PASS BEFORE PASS 1 ********                                                                                                 

		if avg eq 0 then begin   ; start with a monthly-averaged PSF description 
                npre_free=5+ij               ; wcof[0,1], z, offset, amp[0]                                                                                  
                chunk.free_ind[0:4+ij,zp] = 1 ; all free                                                                                                     
        endif
        if avg eq 1 then begin
                npre_free=4+ij               ; wcof[0,1], z, offset, amp[0]                                                                                  
                chunk.free_ind[0:3+ij,zp] = 1 ; all free                                                                                                     
                endif

                ;;;;;;;;;;                                                                                                                                   
        refine_ch = dop_pre(chunk, npre_free, dopenv, lm=lm, demo=demo, pfilt=pfilt, avg=avg)
        chunk=refine_ch
        ;;;;;;;;;;                                                                                                                                           

        ;print,'Zeroth pass, PSF 0: ',(*chunk.free_par)[zp].amp[0]                                                                                           
        if keyword_set(verbose) then begin
                print,'Parameters after Pre-fit, before Pass 1: '
                print,(*chunk.free_par)[zp].wcof[*]
                print,(*chunk.free_par)[zp].z
                print,(*chunk.free_par)[zp].amp[0]
        endif
 
    	chunk.free_ind[0:n_elements(dopenv.psfpix)-1+4+ij,fp] = 1 ; all free
    	chunk.free_ind[4]=0
    	if dopenv.iss_nm eq 'iod' then iod_flag = 1 else iod_flag = 0
    	parinfo=dop_parinfo(chunk, dopenv, pass=pass, iod_flag=iod_flag,osamp=dopenv.osamp)
    	chunk.fit[fp] = 100.
    	functargs = {chunk:chunk, dopenv:dopenv, pass:pass}
    	if keyword_set(verbose) then begin
       		print,'Parameters going to Pass 1 (match above??): '
        	print,parinfo.value
    	endif
	endif  ;**** END SETUP FOR PASS 1 ****

      if pass eq 2 then begin   ;**** BEGIN PASS 2 ****
          if dopenv.n_wcof eq 2 then ifree=[0,1,2,3]
          if dopenv.n_wcof eq 3 then ifree=[0,1,2,3,4]
          iod_flag = 0
          if dopenv.iss_nm eq 'iod' then begin
             iod_flag = 1
             if dopenv.n_wcof eq 2 then ifree=[0,1,3]
             if dopenv.n_wcof eq 3 then ifree=[0,1,2,4]
          endif
          chunk.free_ind[ifree,lp] = 1 ; frz PSF
          parinfo=dop_parinfo(chunk, dopenv, pass=pass, iod_flag=iod_flag, osamp=dopenv.osamp)
          chunk.fit[lp] = 100.
          functargs = {chunk:chunk, dopenv:dopenv, pass:pass, $
                      chunk_arr:chunk_arr}
          if keyword_set(verbose) then begin
             print,'Parameters out of Pass 1: '
             print,(*chunk.free_par)[fp].wcof[*]
             print,(*chunk.free_par)[fp].z
             print,(*chunk.free_par)[fp].amp[*]
             print,'Parameters setup for parinfo (match above??): '
             print,parinfo.value
          endif

       endif  ;**** END SETUP FOR PASS 2 ****

       x = findgen(dopenv.n_pix)    
       y = (*chunk.sobs)
       xfnd=where(chunk_arr.ordt eq chunk.ordt and chunk_arr.pixt eq chunk.pixt,nxfnd) 


    ; FIRST NLLS FITTING - TURN THE (TEMPLATE*IODINE) INTO THE OBSERVATION
    	chi=99.  ;assume the worst 
    	chithresh = 2.5d

		newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, $
			functargs=functargs, maxiter=200, $
			errmsg=errmsg, /iterstop, $
			bestnorm=bestnorm, perror=perror, $
			yfit=syn_fit, weight=wt, status=status,/quiet)

    ; STATUS contains the error_code from dop_fit.pro 
    ;if status lt 0 and status gt -7 then stop

		xnan=where(finite(syn_fit,/nan) ne 0,nxnan)
		
		; NEWPAR LOOKS OK - CHECK WEIGHTS WITH CHAUVENET TEST
		; UPDATE PARINFO.VALUE IF CHISQ FIT IS LT CHITHRESH  
        if n_elements(newpar) eq n_elements(parinfo.value) and $
           nxnan eq 0 then begin
           		npar0=newpar  ; hang on to your hat... 
        		diff = syn_fit - (*chunk.sobs)
        		xpass=chauvenet(diff, npass, reject=failed, nrejects=nf)
       			; REJECT OUTLIERS AND ITERATE WITH CHAUVENET'S CRITERION
       			; DO NOT REJECT MORE THAN 10% OF THE POINTS
				if n_elements(failed) lt 0.1*n_elements(wt) then begin
					wt[failed]=0.d
					newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, $
						functargs=functargs, maxiter=200, /nan, $
						errmsg=errmsg, /iterstop, $
						bestnorm=bestnorm, perror=perror, $
						yfit=syn_fit, weight=wt, status=status,/quiet)
					; unlikely to enter next block
					if n_elements(newpar) lt n_elements(parinfo.value) then begin
						print,'WHAT??  failed after applying chauvenet weights' 
						chi = 98.
						newpar=npar0
					endif
				endif ;chauvenet re-weight	
					
				; ALL IS GOOD, CALCULATE NEW CHISQ AND UPDATE PARINFO	
        		if n_elements(newpar) eq n_elements(parinfo.value) and $
        		   nxnan eq 0 then begin
        				nbad = n_elements(where(wt eq 0.))
						nfree = n_elements(where(parinfo.fixed eq 0))
						dof=n_elements((*chunk.sobs)) - nbad - nfree ; WEIGHT EQ 0 OR FIXED
			  			red_chi_sq = total(  wt * (syn_fit - (*chunk.sobs))^2 )
			  			red_chi_sq = red_chi_sq / dof
			  			chi = sqrt(red_chi_sq) 
			  			if chi le chithresh then parinfo.value = newpar  ; DONE!
			  	endif
		endif   ;chauvenet test

	   ; NEWPAR DOES NOT LOOK OK - ROLL BACK
	   ; OR FIT NOT GREAT, PERTURB THE INITIAL GUESSES
       		if chi gt chithresh and pass eq 1 and perturb eq 1 then begin
       			previous_psfpar1=(*chunk_arr[xfnd-1].free_par)[1].amp   ;17 psfpars from previous chunk
       			previous_psfpar2=(*chunk_arr[xfnd-2].free_par)[1].amp   ;17 psfpars from previous chunk
       			zm1=(*chunk_arr[xfnd-1].free_par)[1].z           ;Doppler shift from previous chunk
       			zm2=(*chunk_arr[xfnd-2].free_par)[1].z           ;Doppler shift from previous chunk
       			zm3=(*chunk_arr[xfnd-3].free_par)[1].z           ;Doppler shift from previous chunk
       			previous_z=double(mean([zm1,zm2,zm3]))                   ;Doppler shift from previous chunk
;print,zm1,zm2,zm3
				initwav = parinfo[0].value
				tries = 0L
				tries1 = 0L
				wavperturb = 0.01   ;2.*parinfo[0].step 
				chibest = 99d
				savparinit = parinfo.value
				savparstep = parinfo.step 
				newparinit = parinfo.value
				newsynfit = y
				newwt=wt
				nbad = n_elements(where(wt eq 0.))
				nfree = n_elements(where(parinfo.fixed eq 0))
				dof=n_elements((*chunk.sobs)) - nbad - nfree ; WEIGHT EQ 0 OR FIXED
				newdof=dof
				endloop = 4d1
				
				; TRY TO FIND A BETTER FIT
				repeat begin
					parinfo[4:20].value=previous_psfpar1
					parinfo[2].value=previous_z
					newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, $
                          	functargs=functargs, maxiter=200, /nan, $
                          	errmsg=errmsg, /iterstop, $
                          	bestnorm=bestnorm, perror=perror, $
                          	yfit=syn_fit, weight=wt, status=status,/quiet)

					    ; loop for catastrophic failure 
          	  		    if n_elements(newpar) eq 1 then begin ; try, try again
             			  	repeat begin   ;27Jun2012 DF
             					parinfo[0].value=initwav + randomn(seed)*wavperturb
								parinfo[4:20].value=previous_psfpar2
             					newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, $
                          				functargs=functargs, maxiter=200, /nan, $
                          				errmsg=errmsg, /iterstop, $
                          				bestnorm=bestnorm, perror=perror, $
                          				yfit=syn_fit, weight=wt, status=quiet)
                          		tries1++
             			  	endrep until n_elements(newpar) eq n_elements(parinfo.value) or tries1 eq endloop 
             			  	if tries1 ge 40 then stop,'screwed!'
             			  	;if tries1 lt 40 then parinfo.value = newpar
          	  		    endif	; catastrophic failure
          	  		    	
          	  		; CALCULATE CHISQ    		  
					nbad = n_elements(where(wt eq 0.))
					nfree = n_elements(where(parinfo.fixed eq 0))
					dof=n_elements((*chunk.sobs)) - nbad - nfree ; WEIGHT EQ 0 OR FIXED
			  		red_chi_sq = total(  wt * (syn_fit - (*chunk.sobs))^2 )
			  		red_chi_sq = red_chi_sq / dof
			  		chi = sqrt(red_chi_sq) 

			  		if chi lt chibest then begin
				 		;save the best solution out of ENDLOOP trials
				 		newparinit = newpar
				 		newsynfit = syn_fit
				 		newwt = wt
				 		newdof = dof	
				 		chibest = chi
			  		endif

			  		if chi ge chithresh then begin 			  		
				 		parinfo[0].value = initwav + randomn(seed)*wavperturb
				 		if (tries mod 8 eq 7) then wavperturb += 0.005  ; bump it up
				 		if (tries mod 15 eq 14) then begin  ;adjust the psf
				 			xmatch=where(chunk_arr.ordt eq chunk.ordt and chunk_arr.pixt eq chunk.pixt,nxmatch)
				 			; grab the psf from the previous chunk
				 			for iq = 4, 20 do newparinit[iq] = (*chunk_arr[xmatch-1].free_par)[1].amp[iq-4]
				 		endif
				 		if (tries mod 20 eq 19) then begin  ;adjust the psf
				 			xmatch=where(chunk_arr.ordt eq chunk.ordt and chunk_arr.pixt eq chunk.pixt,nxmatch)
				 			; grab the wavelength soln from the previous chunk
				 			newparinit[0] = (*chunk_arr[xmatch-1].free_par)[1].wcof[0] + $
				 				dopenv.n_pix*(*chunk_arr[xmatch-1].free_par)[1].wcof[1]
				 		endif
			  		endif
			  		tries++
				endrep until (chi lt chithresh OR tries gt endloop)
				
				newpar = newparinit
				syn_fit = newsynfit
				wt = newwt
				dof=newdof
				red_chi_sq = total(  wt * (syn_fit - (*chunk.sobs))^2 )
				red_chi_sq = red_chi_sq / dof
				chi_iter = sqrt(red_chi_sq) 

				; NEWPAR LOOKS OK - CHECK WEIGHTS WITH CHAUVENET TEST
				; UPDATE PARINFO.VALUE IF CHISQ FIT IS LT CHITHRESH  
        		if n_elements(newpar) eq n_elements(parinfo.value) and $
           			n_elements(syn_fit) eq dopenv.n_pix then begin
           				npar0 = newpar  ; hang on to your hat... 
        				diff = syn_fit - (*chunk.sobs)
        				xpass = chauvenet(diff, npass, reject=failed, nrejects=nf)
       					; REJECT OUTLIERS AND ITERATE WITH CHAUVENET'S CRITERION
       					; DO NOT REJECT MORE THAN 10% OF THE POINTS
						if n_elements(failed) lt 0.1*n_elements(wt) then begin
							wt[failed]=0.d
							newpar=mpfitfun('dop_fit', x, y, parinfo=parinfo, $
								functargs=functargs, maxiter=200, /nan, $
								errmsg=errmsg, /iterstop, $
								bestnorm=bestnorm, perror=perror, $
								yfit=syn_fit, weight=wt, status=status,/quiet)
							; unlikely to enter next block
							if n_elements(newpar) lt n_elements(parinfo.value) then begin
								print,'Pt 2: failed after applying chauvenet weights' 
								newpar=npar0
							endif
						endif ;chauvenet re-weight	
				endif   ;chauvenet test
		   endif ; perturb eq 1		
		
	if n_elements(syn_fit) eq 0 then stop
		
       parinfo.value = newpar
       nbad = n_elements(where(wt eq 0.))
       nfree = n_elements(where(parinfo.fixed eq 0))
       dof=n_elements((*chunk.sobs)) - nbad - nfree ; WEIGHT EQ 0 OR FIXED
       red_chi_sq = total(  wt * (syn_fit - (*chunk.sobs))^2 )
       red_chi_sq = red_chi_sq / dof
       chi = sqrt(red_chi_sq) 
if chi eq 0.0 then chi = 99

       if pass eq 1 and n_elements(newpar) gt 5 then begin
          if dopenv.psfmod eq 'gaussian' then begin
			ip=dop_psf(newpar[4:n_elements(dopenv.psfpix)-1+4], dopenv=dopenv)
          endif;gaussian
          if dopenv.psfmod eq 'bspline' then begin
			ip=dop_psf_bspline(newpar[5:*], dopenv=dopenv, cntr=dopenv.psfcntr)
          endif;bspline
          mod_chunk = dop_chunk_update(chunk, parinfo, ip, syn_fit=syn_fit, chi=chi,$
          	 pass=pass, iod_flag=iod_flag) 
       endif

       if pass eq 2 then begin
          ip=dop_psf_smooth(chunk, chunk.ordob, chunk.pixob, $
                            dopenv=dopenv, chunk_arr=chunk_arr)
          mod_chunk = dop_chunk_update(chunk, parinfo, ip, syn_fit=syn_fit, chi=chi,$ 
          	pass=pass, iod_flag=iod_flag)
       endif

       if keyword_set(verbose) then begin
          print,'After Pass 1 chunk free_par values (match above???): '
          print,(*mod_chunk.free_par)[pass].wcof[*]
          print,(*mod_chunk.free_par)[pass].z
          print,(*mod_chunk.free_par)[pass].amp[*]
       endif

       if ~keyword_set(demo) then $
          print,'Pass: ', strcompress(string(pass),/remove_all), $
                '  Chunk: ',strcompress(string(chunk.ordob),/remove_all), '  ',$
                strcompress(string(chunk.pixob),/remove_all),'  Red chi: ', $
                strmid(strcompress(string(chi),/remove_all),0,4)
    
    	if chi eq 0.0 then begin  ; this is happening in the tries++ loop
    		print, wt
			;print, (*mod_chunk.free_par)[1].wcof[0]-winit, '  ',(*mod_chunk.free_par)[1].wcof[1]-dinit
;			stop
		endif
		;if chi gt 5 then stop
		
       if keyword_set(demo) then begin 
          print,'Pass: ', strcompress(string(pass),/remove_all), $
                '  Chunk: ',strcompress(string(chunk.ordob),/remove_all), '  ',$
                strcompress(string(chunk.pixob),/remove_all),'  Red chi: ', $
                strmid(strcompress(string(chi),/remove_all),0,4)

        ; FIRST PANEL: DECONV SPECTRUM AND FTS IODINE
          !p.charsize=1.6
          !x.charsize=1.8
          !y.charsize=1.8
          !p.multi=[0,1,3]
          !x.omargin=[8,2]
          !y.omargin=[2,2]
          if dopenv.osamp eq 4 then begin
             cent = 60
             yrr=[0,0.8]
          endif
          if dopenv.osamp eq 8 then begin
             cent=56
             yrr=[0,2]
          endif

          wav_coef = (*mod_chunk.free_par)[pass].wcof
          wobs = poly(findgen(dopenv.n_pix), wav_coef)
          xr1=max( [min(*mod_chunk.wis), min(wobs)] )
          xr2=min( [max(*mod_chunk.wis), max(wobs)] )
          plot, ip, col=1, title='IP', /xsty,xra=[40,80], /ysty, yra=yrr
          plots, [cent,cent], yrr
          plot, (*mod_chunk.wis), (*mod_chunk.sis),col=1,title='!6 ISS and I2', $
                xra=[xr1,xr2],/xsty, yra=[0,1.1],/ysty, yticks=2,ytickv=[0, 0.5, 1.0],$
                ytickname=['0','0.5','1.0']
          oplot, (*mod_chunk.wiod), (*mod_chunk.siod), col=200


       ; SECOND PANEL: MODEL AND OBSERVED SPECTRUM
          plot, wobs, (*mod_chunk.sobs), xra=[xr1,xr2], /xsty, yticks=3,$
               /ysty, xtit='!6 Wavelength', $
               titl='!6 Model (red) and Obs Spectrum'
          !p.color=222
          oplot, wobs, syn_fit
          if nf gt 0 then oplot, wobs[failed], syn_fit[failed], col=130, ps=7, symsize=1.4,thick=2
          oploterr, wobs, syn_fit, sig

          !p.color = 1
       endif  ; demo

       return, mod_chunk 

end

