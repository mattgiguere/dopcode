pro dr_iodkeck_loop, tag=tag, run=run, demo=demo, verbose=verbose, $
                     jdmin=jdmin, jdmax=jdmax, $
                     observatory=observatory, obmin=obmin,obmax=obmax

if ~keyword_set(tag) then tag='s'
observatory='keck'
runlen=strlen(run)   ; eg run='rj101'
deck_select='C2'

; e.g., dr_iod_loop,tag='b',observatory='ctio',/demo

  if ~keyword_set(obmin) then obmin=0
  if ~keyword_set(obmax) then obmax=9999

  if observatory eq 'keck' then begin
     restore,'/tous/mir3/bary/kbcvel.dat'
     xmatch=where(strmid(bcat.obsnm,0,runlen) eq run and $
                  strmid(bcat.objnm,0,2) eq 'HR' and $                  
                  bcat.obtype eq 'o',num)  ;gets all deckers 
  endif

   if num eq 0 then stop 
   obs_id=bcat[xmatch].obsnm
   flag=intarr(n_elements(obs_id))
   for i=0,n_elements(obs_id)-1 do begin
      hd=headfits('/tous/mir3/fitspec/'+obs_id[i]+'.fits')
      decker=sxpar(hd,'DECKNAME')
      if strcompress(decker,/rem) eq deck_select then flag[i]=1 else flag[i]=0
   endfor

obs_id=obs_id[where(flag eq 1)]
num=n_elements(obs_id) 

   for jobs = 0, num-1 do begin 
      if observatory eq 'keck' then dopenv=keck_init(obs_id[jobs], 'iod', 0.0, iss_obnm='iod',tag=tag) 
         fck=findfile(dopenv.obs_dir+dopenv.obsnm+'.fits',count=count)
      if count eq 0 then print,'FITS FILE NOT FOUND: ',dopenv.obsnm 
      if count eq 0 then goto, JUMP
;         restore,dopenv.bary_file  ;bcat
;         x=where(bcat.obsnm eq obs_id[jobs],nx) 
;      if nx gt 0 then begin
;         bcat=bcat[x[0]]
;         objnm = bcat.objnm
         found=findfile(dopenv.files_dir+'cd'+tag+'b*'+obs_id[jobs], count=runcount)
         if runcount eq 1 then print, found
         if runcount eq 0 then begin
            dop_main, dopenv.obj_nm, dopenv, obsnm=obs_id[jobs], $
                              pass=1, tag=tag, /iod_soln, lm='mpf', demo=demo 
            dop_main, dopenv.obj_nm, dopenv, obsnm=obs_id[jobs], pass=2, lm='mpf', tag=tag, /iod_soln 
            dop_create_ipcf,tag=tag,observatory=observatory,goodfit=5.
            print, 'Done: '+found
         endif
;      endif
      JUMP: 
      if count eq 0 then print,'Iodine obs not found to run'
   endfor   ;j loop of obs_id

end
