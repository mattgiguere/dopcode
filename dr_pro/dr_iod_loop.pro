pro dr_iod_loop, tag=tag, run=run, demo=demo, verbose=verbose, jdmin=jdmin, jdmax=jdmax, $
                 observatory=observatory, obmin=obmin,obmax=obmax, obarr=obarr 

;fischer runs the iodine modeling 
; added obsnm, which is an array 
; e.g., dr_iod_loop,tag='rchi',run='rchi111102',observatory='ctio4k',obarr=['1245', '1246', '1247']
; dr_iod_loop,tag='dfs',run='ri67',observ='lick'
; dr_iod_loop,tag='dfg',run='rg92',observ='lick', obmin=216, obmax=233 

if ~keyword_set(observatory) then begin
   observatory=''
   read,'enter observatory (ctio, ctio4k, lick, keck): ', observatory
endif

if ~keyword_set(tag) then begin
   tag=''   
   read,'enter observatory (e.g., b): ', tag
endif

;if ~keyword_set(run) then begin
;   run=''
;   read,'enter observatory (e.g., rqa13.1200): ', run
;endif

; e.g., dr_iod_loop,tag='b',observatory='ctio',/demo

  if ~keyword_set(obmin) then obmin=0
  if ~keyword_set(obmax) then obmax=9999

  if observatory eq 'ctio' then begin
     restore,'/mir7/bary/qbcvel.dat'
     xmatch=where(strmid(bcat.obsnm,0,5) eq run and $
                  strmid(bcat.objnm,0,2) eq 'io' and bcat.obtype eq 'i' and $;bcat.jd gt jdmin,num)
                  fix(strmid(bcat.obsnm,6,4)) ge obmin and fix(strmid(bcat.obsnm,6,4)) le obmax,num)
  endif

  if observatory eq 'ctio4k' and keyword_set(obmin) then begin
     restore,'/mir7/bary/qbcvel.dat'
;     xmatch=where(strmid(bcat.obsnm,0,5) eq run and $
;                  strmid(bcat.objnm,0,2) eq 'HR' and bcat.obtype eq 'o' and $;bcat.jd gt jdmin,num)
;                  fix(strmid(bcat.obsnm,6,4)) ge obmin and fix(strmid(bcat.obsnm,6,4)) le obmax,num)
     xmatch=where(strmid(bcat.obsnm,0,10) eq run and $ ;strmid(bcat.objnm,0,2) eq 'HR' and 
                       bcat.obtype eq 'i' and $;bcat.jd gt jdmin,num)
                  fix(strmid(bcat.obsnm,11,4)) ge obmin and fix(strmid(bcat.obsnm,11,4)) le obmax,num)
  endif

  if observatory eq 'ctio4k' and keyword_set(obarr) then begin
     restore,'/mir7/bary/qbcvel.dat'
     xmatch=lonarr(n_elements(obarr))
     for i=0,n_elements(obarr)-1 do begin 
     	xm=where(strmid(bcat.obsnm,0,10) eq run and bcat.obtype eq 'i' and $
                  fix(strmid(bcat.obsnm,11,4)) eq obarr[i], nxm)
        if nxm gt 0 then xmatch[i]=xm
     endfor
     s=where(xmatch gt 0,nx) 
     xmatch=xmatch[s]
     num=nx
  endif

  if observatory eq 'keck' then begin
     runlen=strlen(run)
     restore,'/tous/mir3/bary/kbcvel.dat'
     xmatch=where(strmid(bcat.obsnm,0,runlen) eq run and strmid(bcat.objnm,0,2) eq 'HR' and bcat.obtype eq 'o',num)
  endif

  if observatory eq 'lick' then begin
     restore,'/mir1/bary/bcvel.dat'
     if keyword_set(run) then xmatch=where(strmid(bcat.obsnm,0,4) eq run and $
     	strupcase(strmid(bcat.objnm,0,2)) eq 'HR' and bcat.obtype eq 'o' and strmid(bcat.obsnm,5,3) ge obmin $
     	and strmid(bcat.obsnm,5,3) le obmax,num) else $ ;and bcat.jd gt 15434.5
     xmatch=where(bcat.jd gt jdmin and bcat.jd lt jdmax and strmid(bcat.objnm,0,2) eq 'HR' and bcat.obtype eq 'o',num)
  endif

   if num eq 0 then stop 
   obs_id=bcat[xmatch].obsnm

   for jobs = 0, num-1 do begin 
      if observatory eq 'ctio' then dopenv=ctio_init(obs_id[jobs], 'iod', 0.0, iss_obnm='iod',tag=tag) 
      if observatory eq 'ctio4k' then dopenv=ctio4k_init(obs_id[jobs], 'iod', 0.0, iss_obnm='iod',tag=tag) 
      if observatory eq 'keck' then dopenv=keck_init(obs_id[jobs], 'iod', 0.0, iss_obnm='iod',tag=tag) 
      if observatory eq 'lick' then dopenv=lick_init(obs_id[jobs], 'iod', 0.0, iss_obnm='iod',tag=tag) 
         fck=findfile(dopenv.obs_dir+dopenv.obsnm+'.fits',count=count)
      if count eq 0 then print,'FILE NOT FOUND: ',dopenv.obsnm 
      if count eq 0 then goto, JUMP
         restore,dopenv.bary_file  ;bcat
         x=where(bcat.obsnm eq obs_id[jobs],nx) 
      if nx gt 0 then begin
         bcat=bcat[x[0]]
         objnm = bcat.objnm
         found=findfile(dopenv.files_dir+'cd'+tag+'b*'+obs_id[jobs], count=count)
         if count eq 1 then print, found
         if count eq 0 then begin
            dop_main, dopenv.obj_nm, dopenv, obsnm=obs_id[jobs], $
                              pass=1, tag=tag, /iod_soln, lm='mpf', demo=demo 
            dop_main, dopenv.obj_nm, dopenv, obsnm=obs_id[jobs], pass=2, lm='mpf', tag=tag, /iod_soln 
            dop_create_ipcf,tag=tag,observatory=observatory,goodfit=2.
            print, 'Done: '+found
         endif
      endif
      JUMP: 
      if nx eq 0 then print,'Iodine obs not found to run'
   endfor   ;j loop of obs_id

end
