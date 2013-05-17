PRO DOP_VEL, STARNAME, LABEL, MAXCHI=MAXCHI, MINCTS=MINCTS

; Read all vd structures and weight the pixels
; Ala vdcube, vdclean (by GWM, RPB) 
; 
; INPUT
; starname: 'HIP8102'
; label: 'b'
; maxchi: 3.0 (default is 2.0)
; 
; OUTPUT
; cube_array: an array with weighted chunk information
; 
; Fischer, SFSU, Apr 2008

  loadct,39
  !p.background=255
  !p.color=1

   vstnm='vstbank/vst'+strlowcase(starname)+'.dat'
   c_light = 2.99792458d8
   restore,'/mir1/bary/bcvel.dat'
   fname=strcompress('vd'+label+starname+'*', /remove_all)
   files=findfile('/mir1/files/'+fname, count=nobs)
   flag=intarr(nobs)
   restore,files[0]     ; vd
   nchunk=n_elements(vd.ordt)

   ; SETUP THE CUBE AND STRUCTURES
   restore,files[0]       ; first one found
   st=strpos(files[0],'_')
   obnm=strmid(files[0], st+1, strlen(files[0])-1)

   if ~keyword_set(maxchi) then maxchi=2.0
   if ~keyword_set(mincts) then mincts=1000.
   if ~keyword_set(maxcts) then maxcts=100000.

   form1 = '(a16, f12.2, i15, f10.2)'
   print,' '
   print,'----------------------------------------------------------------'
   print,' Restored File  Median Velocity   Median Photons  Median'
   print,'                 *All Chunks*      *All Chunks*    Fit'
   print,'----------------------------------------------------------------'

   st=strpos(files[0],'_')
   obsnm=strmid(files[0], st+1, strlen(files[0])-1)

 ; DEFINE THE CUBE STRUCTURE TAG NAMES
   cube_one = {obnm:obnm,                 $ ;observation name
               z:0.0d,                    $ ;doppler z = v/c = dl/l
               bc:0.0,                    $ ;barycentric correction
               vel:fltarr(nchunk),        $ 
               fit:fltarr(nchunk),        $
               weight:fltarr(nchunk),     $ ;initially, the IS weight
               cts:fltarr(nchunk)}
   cube=replicate(cube_one, nobs)

   gdind=[0]  ;track the good observations
;   velarr=fltarr(nobs)

;LOOP THROUGH ALL OBSERVATIONS - REJECT BAD OBSERVATIONS
   for i=0, nobs-1 do begin
;      print, 'Weighting the observations'
      restore,files[i]
      st=strpos(files[i],'_')
      obsnm=strmid(files[i], st+1, strlen(files[i])-1)
      cube[i].obnm=obsnm
      xx=where(obsnm eq bcat.obsnm,nxx)
      if nxx ne 1 then stop
      cube[i].bc = bcat[xx].bc

    ; FILL THE CUBE STRUCTURE
      for k=0,nchunk-1 do begin
         cube[i].z = double(vd[k].z)
       ; ADD IN THE BC VELOCITY
         cube[i].vel[k]=double(vd[k].z*c_light + cube[i].bc)
         cube[i].fit[k]=vd[k].fit
         cube[i].weight[k] = vd[k].weight ; from IS
;         if cube[i].weight[k] le 0.0 then cube[i].weight[k]=0.0
         cube[i].cts[k] = vd[k].cts
      endfor ;k
      medvel = median(cube[i].vel) ;second pass vel
      medcts = long(median(cube[i].cts))
      medchi = median(cube[i].fit) ;second pass chi fit

    ; IDENTIFY BAD OBSERVATIONS
      message_1=''   &   message_2=''
      if medchi gt maxchi then message_1 = 'High Chi-sq: '
      if medchi eq 0.0 then message_1 = 'Chi-sq = 0.0 '
      if medcts lt mincts then message_2 = 'Low Counts: ' 
      if medcts gt maxcts then message_2 = 'High Counts: '
      if message_1 ne '' or message_2 ne '' then begin
         cube[i].weight=0.0
         print,message_1, message_2
      endif
      if message_1 eq '' and message_2 eq '' then gdind=[gdind,i]
   endfor   ; filling cube 

 ; REJECT BAD OBSERVATIONS - RESIZE THE CUBE ARRAY AND UPDATE NOBS
   ngd = n_elements(gdind)-1  ;trim off the dummy start
   gdind = gdind[1:ngd]
   cube = cube[gdind,*]
   nobs = ngd

 ; LOOP THROUGH ALL CHUNKS - REJECT BAD CHUNKS
   print, 'Weighting the chunks'
   err=fltarr(nobs,nchunk)
   fit=fltarr(nobs,nchunk)
   vel=fltarr(nobs,nchunk)
   rvel=fltarr(nobs,nchunk)
   
   sig=fltarr(nchunk)

   for i=0,nobs-1 do begin
      for k=0,nchunk-1 do begin
         if double(cube[i].weight[k]) eq 0.0 then err[i,k]=99. else $
            err[i,k]=1./sqrt(cube[i].weight[k])     
         fit[i,k]=cube[i].fit[k]
         vel[i,k]=cube[i].vel[k]
      endfor ;k
      mvel = median(cube[i].vel[*])
      rvel[i,*] = cube[i].vel[*] - mvel 
      gd_ind=chauvenet(rvel[i,*],npass,reject=failed,/iterate)
      cube[i].weight[failed]=0.0
      err[i,failed]=99.
;      if i eq 0 then plot,rvel[i,gd_ind],col=1,ps=8,symsize=0.6, $
;                          title='Chauvenet selected RVs' else $
;               oplot,rvel[i,gd_ind],col=i*2.2, symsize=0.6,ps=8
   endfor

;plothist,cube.weight,bin=0.001,title='Distrib of chunk weights',xtit='Weights'
;wait,2

   mederr = median(err)         ;second pass photon weighting
   medfit = median(fit)         ;second pass chi fit


 ; DISCARD POOR CHI-SQ FITS AND LARGE PHOTON ERRORS
   ierr=sort(err) & s_err=err(ierr)
   max_err=s_err[0.997*nchunk*nobs] < 99.   ; 3-sigma
   max_err=s_err[0.954*nchunk*nobs] < 99.   ; 2-sigma
   ifit=sort(fit) & s_fit=fit(ifit)
   max_fit=s_fit[0.997*nchunk*nobs]         ; 3-sigma
   max_fit=s_fit[0.954*nchunk*nobs]         ; 2-sigma

   if keyword_set(maxchi) then max_fit=min(maxchi, max_fit)

print,max_err
print,max_fit
plothist,err,bin=20,xra=[0,200],title='Distrib of Err',xtit='Chunk err [m/s]'
wait,2
plothist,fit,bin=0.2,xra=[0,5],title='Distrib of Fits',xtit='Chunk chi fits'
wait,2

 ; SET BAD CHUNK WEIGHTS TO ZERO
   for i=0,nobs - 1 do begin
      for k=0,nchunk-1 do begin
         if err[i,k] ge max_err then cube[i].weight[k]=0.0
         if fit[i,k] ge max_fit then cube[i].weight[k]=0.0
      endfor    
   endfor
;plothist,cube.weight,bin=0.001,title='Distrib Weights',xtitle='Weights'
;wait,2
; ; WITHIN A CHUNK SET, WEIGHT SETS ACCORDING TO VELOCITY RMS
   fudge=1.0
   for k=0,nchunk-1 do begin
      igd = where(double(cube[*].weight[k]) gt 0.0d,ngd) ; igd = fltarr(nobs)
      if ngd gt 3 then sig[k] = stddev(rvel[igd,*]) else sig[k]=1000.   
   endfor
   for i=0,nobs-1 do begin
    ; RATIO EACH CHUNK RMS TO THE MEAN RMS
      const = median(abs(rvel[i,igd])/sig[igd])
      sig_obs = const*sig*fudge
      for k=0,nchunk-1 do cube[i].weight[k] = 1./sig_obs[k]^2.
   endfor
;stop
 ; DEFINE THE CF STRUCTURE TAG NAMES
   cf = {obnm:obnm,   $         ;observation name
         bc:0.0,           $    ;barycentric correction
         jd:0.0d,          $    ;jd of observation
         dewar:24,         $    ;dewar number
         gain:6.1,         $    ;dewar gain
         cts:0L,           $    ;photon counts
         mnvel:0.0,        $    ;mean vel
         mdvel:0.0,        $    ;median vel
         errvel:0.0,       $    ;doppler error
         z:0.0d,           $    ;doppler z
         mdchi:0.0,        $    ;median chi fit
         nchunk:nchunk}         ;number of chunks

   cf=replicate(cf, nobs)

 ; FILL THE CF STRUCTURE
   for ii=0,nobs - 1 do begin
      cf[ii].obnm=cube[ii].obnm
      xb=where(cf[ii].obnm eq bcat.obsnm,nxb)
      if nxb ne 1 then stop else begin
         cf[ii].jd=bcat[xb].jd
         cf[ii].bc = bcat[xb].bc 
      endelse
      velobs = reform(cube[ii].vel)  ; velobs=fltarr(nchunk)
    ; INDENTIFY INDICES OF CHUNKS WITH GOOD WEIGHTS AND LOW SCATTER
      resid = abs(velobs - median(velobs))
      ires=sort(resid)
      num_res=n_elements(resid)
      thresh = 0.997   ; 3 sigma
      thresh = 0.954   ; 2 sigma
      gd1 = where(cube[ii].weight lt 100 and double(cube[ii].weight) gt 0.0)
      gd2 = ires(0:thresh*(num_res-1))
      gd = [gd1, gd2]
      m=sort(gd)  &  gd=gd[m]
      g=uniq(gd)  &  gd=gd[g]

      velobs = cube[ii].vel[gd]
      wt = cube[ii].weight[gd]
      fit = cube[ii].fit[gd]
      cf[ii].mdchi = median(fit)                 ; median fit
      cf[ii].cts = median(cube[ii].cts[gd])     ; median counts
      cf[ii].mdvel = median(velobs)              ; median vel
      cf[ii].mnvel = total(velobs*wt)/total(wt)  ; weighted mean vel
      cf[ii].mnvel = cf[ii].mnvel - cf[ii].z*cf[ii].bc  ; new correction
      cf[ii].errvel = 1./sqrt(total(wt))
   endfor

cf1=cf
cf3=cf
save,cf1,cf3,f=vstnm
if keyword_set(patch) then patch_d8, star=strlowcase(starname)
restore,'vstbank/vst'+strlowcase(starname)+'.dat'
velplot,cf3

stop
end
