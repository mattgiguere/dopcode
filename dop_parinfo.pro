FUNCTION DOP_PARINFO, chunk, dopenv, pass=pass, iod_flag=iod_flag, osamp=osamp, avg=avg

 ; CREATE A TEMPORARY STRUCTURE FOR FREE PARAMETERS
  
c_light= dopenv.c_light

;;;;;;;;;;;;;;;;PASS = 1;;;;;;;;;;;;;;;;;;;;;;;;;;;
if pass eq 1 then begin 

if dopenv.n_wcof eq 3 then names=['WCOF0','WCOF1','WCOF2','Z','OFFSET','AMP0','AMP1','AMP2','AMP3',$
          'AMP4','AMP5','AMP6','AMP7','AMP8','AMP9','AMP10','AMP11','AMP12',$
          'AMP13','AMP14','AMP15','AMP16']  else $
          
          names=['WCOF0','WCOF1','Z','OFFSET','AMP0','AMP1','AMP2','AMP3',$
          'AMP4','AMP5','AMP6','AMP7','AMP8','AMP9','AMP10','AMP11','AMP12',$
          'AMP13','AMP14','AMP15','AMP16'] 


 ; START FITTING WITH PAR'S FROM THE PREVIOUS PASS
   par=double( [ (*chunk.free_par)[pass-1].wcof[0],  $
        (*chunk.free_par)[pass-1].wcof[1],  $
        (*chunk.free_par)[pass-1].z*c_light,$
        (*chunk.free_par)[pass-1].offset,   $ 
        (*chunk.free_par)[pass-1].amp[0],   $
        0.001, 0.001, 0.001,                $
;             dopenv.psfsig[0], $
        (*chunk.free_par)[pass-1].amp[4],   $
        (*chunk.free_par)[pass-1].amp[5],   $
        (*chunk.free_par)[pass-1].amp[6],   $
        (*chunk.free_par)[pass-1].amp[7],   $
        (*chunk.free_par)[pass-1].amp[8],   $  ; near center
        (*chunk.free_par)[pass-1].amp[9],   $  ; near center
        (*chunk.free_par)[pass-1].amp[10],  $
        (*chunk.free_par)[pass-1].amp[11],  $
        (*chunk.free_par)[pass-1].amp[12],  $
        (*chunk.free_par)[pass-1].amp[13],  $
        0.001, 0.001, 0.001])
  
   npar=n_elements(names)

 ; SETUP THE PARINFO STRUCTURE
   parinfo = {value: 0.0d,        $    ; double precision  
              fixed: 0,          $
              limited: [0,0],    $ ; use with caution 
              limits: fltarr(2), $
              parname: '?',      $
              step: 0.0d,        $
              relstep: 0.00,     $
              mpside: 2}       ; 0, 1, -1, 2   
   parinfo=replicate(parinfo, npar)

   for i=0,npar-1 do begin
      parinfo[i].parname = names[i]
      parinfo[i].value = par[i]
   endfor

   
  ; SPECIAL PARAMETER INSTRUCTIONS     
  ; PAR[0]: WAVELENGTH                    
    parinfo[0].step = 0.001d
	parinfo[0].limited=[1,1]   ;df 11/13/11
	parinfo[0].limits[0]=parinfo[0].value-0.1
	parinfo[0].limits[1]=parinfo[0].value+0.1
	
  ; PAR[1]: DISPERSION
    parinfo[1].step = 0.0001d
    parinfo[1].limited = [1,1]
    parinfo[1].limits[0] = par[1] - 0.005  ; *(1.-0.01)
    parinfo[1].limits[1] = par[1] + 0.005  ; *1.01
    if dopenv.n_wcof eq 3 then begin
    	ij=1 
    	parinfo[1+ij].step=0.00001
    	parinfo[1+ij].limited=[1,1]
    	parinfo[1+ij].limits=parinfo[1].limits/10.
    endif 
    if dopenv.n_wcof eq 2 then ij=0


  ; PAR[2]: DOPPLER SHIFT        
    parinfo[2+ij].step = 5.d  
    if iod_flag eq 1 then begin 
       parinfo[2+ij].value = 0.0d 
       parinfo[2+ij].fixed = 1
    endif
                     
  ; PAR[3]: SCALE CONTINUUM (THIS USED TO BE SCATTERED LIGHT)  
    parinfo[3+ij].step = 0.01d
    parinfo[3+ij].limited = [1,1]
    parinfo[3+ij].limits = [0.9d, 1.1d]
    
if dopenv.psfmod eq 'gaussian' then begin   
  ; PAR[4]: CENTRAL GAUSSIAN 
 ;   parinfo[4].value=0.9d
    parinfo[4].step=0.001
    parinfo[4+ij].fixed=0          ; width
    
  ;[1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16]
  ;[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
  
  ; PAR[5] - PAR[20]: DEFINE SIDE GAUSSIANS   
      for i=5,n_elements(par)-1 do begin
         parinfo[i+ij].step = 0.001
      endfor

  ; CENTRAL FLANKING GAUSSIANS MUST BE LOWER THAN CENTRAL GAUSSIAN
  ; THIS CRASHES THE DOP ANALYSIS! COMPARES AMPLITUDES TO WIDTH OF
  ; CENTRAL GAUSSIAN
  ;  parinfo[12+ij:13+ij].limited=[0,1]
  ;  parinfo[12+ij:13+ij].limits[1] = 0.95*parinfo[4+ij].value 
     
;      nf=n_elements(parinfo)  ; 4+17 = 21
;       parinfo[5+ij:7+ij].limited=[1,1]
;      parinfo[18+ij:20+ij].limited=[1,1] 
;      parinfo[5+ij].limits = [-0.05,0.05]
;      parinfo[6+ij].limits = [-0.09,0.09]
;      parinfo[7+ij].limits = [-0.1,0.14]
;      parinfo[18+ij].limits = [-0.1,0.14]
;      parinfo[19+ij].limits = [-0.09,0.09]
;      parinfo[20+ij].limits = [-0.05,0.05]
;;;      parinfo[5:20].limited=[1,1]
;;;      parinfo[5:20].limits=[-0.2,0.8]
 endif  ; psfmod=gaussian

 if dopenv.psfmod eq 'bspline' then begin
 	; THESE PARAMETERS ARE not USED IN THE BSPLINE MODEL
 		parinfo[4].fixed=1
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

      for i=5,n_elements(par)-1 do begin
         parinfo[i].step = 0.001
         parinfo[i].limits = [-0.1,0.8]
      endfor
      
   for i=5,npar-1 do parinfo[i].value = psfpar[i-5]
 endif  ;psfmod=bspline


   endif  ; pass eq 1


;;;;;;;;;;;;;;;;PASS = 2;;;;;;;;;;;;;;;;;;;;;;;;;;;
if pass eq 2 then begin 

   names=['WCOF0','WCOF1','Z','OFFSET']

 ; USE UPDATED PARS FROM PASS 1 

   par=double([(*chunk.free_par)[pass-1].wcof[0], $
        (*chunk.free_par)[pass-1].wcof[1], $  
        (*chunk.free_par)[pass-1].z*c_light, $
        (*chunk.free_par)[pass-1].offset])
   nfree=n_elements(names)

 ; SETUP THE PARINFO STRUCTURE
   parinfo={value: 0.0d,        $    ; double precision  
             fixed: 0,          $
             limited: [0,0],    $    ; use with caution 
             limits: fltarr(2), $
             parname: '?',      $
             step: 0.0d,        $
             relstep: 0.00,     $
             mpside: 2}              ; 0, 1, -1, 2   
   parinfo=replicate(parinfo, nfree)

   for i=0,nfree-1 do begin
      parinfo[i].parname = names[i]
      parinfo[i].value = par[i]
   endfor
   
  ; SPECIAL PARAMETER INSTRUCTIONS     
  ; PAR[0] AND PAR[1]: WAVELENGTH AND DISPERSION                                            
    parinfo[0].step = 0.001d
	parinfo[0].limits[0]=parinfo[0].value-0.3
	parinfo[0].limits[1]=parinfo[0].value+0.3
    parinfo[1].step = 0.0001d
        parinfo[1].limited = [1,1]
        parinfo[1].limits[0] = par[1] - 0.001  ; *(1.-0.01)
        parinfo[1].limits[1] = par[1] + 0.001  ; *1.01   

   if dopenv.n_wcof eq 3 then begin
    	ij=1 
    	parinfo[1+ij].step=0.0001
 ;   	parinfo[1+ij].limited=[1,1]
 ;   	parinfo[1+ij].limits=parinfo[1].limits/10.
   endif 
   if dopenv.n_wcof eq 2 then ij=0

                                
  ; PAR[2]: DOPPLER SHIFT        
   parinfo[2+ij].step = 5.   
   if iod_flag eq 1 then parinfo[2+ij].fixed = 1

                              
  ; PAR[3]: SCALE CONTINUUM (THIS USED TO BE SCATTERED LIGHT)                                      
   parinfo[3+ij].step = 0.01d
   parinfo[3+ij].limited = [1,1]
   parinfo[3+ij].limits = [0.95d, 1.05d]
endif

;;;;;;;;;;;;;;;;PASS = 3;;;;;;;;;;;;;;;;;;;;;;;;;;;
if pass eq 3 then begin 

   names=['WCOF0','WCOF1','Z','OFFSET']

 ; USE PARS FROM PASS 1 

   par=double([(*chunk.free_par)[pass-3].wcof[0], $
        (*chunk.free_par)[pass-3].wcof[1], $  
        (*chunk.free_par)[pass-3].z*c_light, $
        (*chunk.free_par)[pass-3].offset])
   nfree=n_elements(names)

 ; SETUP THE PARINFO STRUCTURE
   parinfo={value: 0.0d,        $    ; double precision  
             fixed: 0,          $
             limited: [0,0],    $    ; use with caution 
             limits: fltarr(2), $
             parname: '?',      $
             step: 0.0d,        $
             relstep: 0.00,     $
             mpside: 1}              ; 0, 1, -1, 2   
   parinfo=replicate(parinfo, nfree)

   for i=0,nfree-1 do begin
      parinfo[i].parname = names[i]
      parinfo[i].value = par[i]
   endfor
   
  ; SPECIAL PARAMETER INSTRUCTIONS     
  ; PAR[0] AND PAR[1]: WAVELENGTH AND DISPERSION                                            
    parinfo[0].step = 0.001d
;   parinfo[0].limited = [1,0]
;   parinfo[0].limits[0] = 0.0

   parinfo[1].step = 0.0001d
   
   if dopenv.n_wcof eq 3 then begin
    	ij=1 
    	parinfo[1+ij].step=0.00001
;    	parinfo[1+ij].limited=[1,1]
;    	parinfo[1+ij].limits=parinfo[1].limits/10.
   endif 
   if dopenv.n_wcof eq 2 then ij=0

                                
  ; PAR[2]: DOPPLER SHIFT        
   parinfo[2+ij].step = 0.01d   
   if iod_flag eq 1 then parinfo[2+ij].fixed = 1

                              
  ; PAR[3]: SCALE CONTINUUM (THIS USED TO BE SCATTERED LIGHT)                                      
   parinfo[3+ij].step = 0.001d
;   parinfo[3+ij].limited = [1,1]
;   parinfo[3+ij].limits = [0.95d, 1.05d]
endif

   return, parinfo

end

