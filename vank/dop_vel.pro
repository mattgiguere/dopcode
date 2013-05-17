PRO DOP_VEL, STARNAME, LABEL, MAXCHI=MAXCHI, MINCTS=MINCTS, $
             C_ARRAY

; Read all vd structures and weight the pixels
; Ala vdcube, vdclean (by G. Marcy) 

; INPUT
; starname: 'HIP8102'
; label: 'b'
; maxchi: 3.0 (default is 2.0)
; 
; OUTPUT
; c_array: an array of chunk_arr that are deemed good
; 
; Fischer, SFSU, Apr 2008

   c_light = 2.99792458d8
   if ~keyword_set(label) then label='b'
   fname=strcompress('vd'+label+starname+'*', /remove_all)
   files=findfile('/mir7/files/'+fname, count=nobs)
   flag=intarr(nobs)
   obsnm=strarr(nobs)
   for i=0,nobs-1 do begin
      st=strpos(files[i],'_')
      obsnm[i]=strmid(files[i], st+1, strlen(files[i])-1)
   endfor

   cf=dop_mkcf(nobs,obsnm)
   if ~keyword_set(maxchi) then maxchi=2.0
   if ~keyword_set(mincts) then mincts=1000.
   if ~keyword_set(maxcts) then maxcts=100000.

   form1 = '(a16, f12.2, i15, f10.2)'

;PRINT HEADER
   print,' '
   print,'----------------------------------------------------------------'
   print,' Restored File  Median Velocity   Median Photons  Median'
   print,'                 *All Chunks*      *All Chunks*    Fit'
   print,'----------------------------------------------------------------'

;LOOP THROUGH ALL OBSERVATIONS
   for i=0, nobs-1 do begin
      restore,files[i]

      nck=n_elements(chunk_arr)  & vel=fltarr(nck)  & fit=fltarr(nck)
      for k=0,nck-1 do begin
         vel[k]=(*chunk_arr[k].free_par)[1].z*c_light + chunk_arr[k].bccor
         fit[k]=(*chunk_arr[k].output)[1].fit
      endfor

      medvel = median(vel) ;second pass vel
      medcts = long(median(chunk_arr.cts))
      medchi = median(fit) ;second pass chi fit
      message_1=''   &   message_2=''
      if medchi gt maxchi then message_1 = 'High Chi-sq: ';,strtrim(medvel,2)
      if medchi eq 0.0 then message_1 = 'Chi-sq = 0.0 '
      if medcts lt mincts then message_2 = 'Low Counts: ';,strtrim(medcts,2)
      if medcts gt maxcts then message_2 = 'High Counts: ';,strtrim(medcts,2)

      if message_1 eq '' and message_2 eq '' then flag[i]=1 else flag[i]=0
      print,obsnm[i], medvel, medcts, medchi, format=form1
   endfor

      
;   xgd=where(flag eq 1,nxgd)
;   if nxgd gt 1 then begin
;      restore,files[xgd[0]]     ;first good chunk_arr
;      c_array = replicate({chunk_arr}, nxgd)
;      for i = 0,nxgd-1 do begin
;         print,files[i], medvel, medcts, medchi, format=form1
;         restore,files[xgd[i]]  ;first good chunk_arr
;         if i eq 0 then c_array[i] = chunk_arr else c_array[i] = [chunk_arr, c_array[i]]
;      endfor
;endif


end
