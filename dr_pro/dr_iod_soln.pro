pro dr_iod_soln, obsnm, tag=tag, demo=demo, verbose=verbose, psfmod=psfmod, $
	ctio=ctio,ctio4k=ctio4k, keck=keck,lick=lick, date=date, yrmo=yrmo, avg=avg, $
	vdavgnm=vdavgnm, cdnear_name=cdnear_name, shft_style=shft_style

; fischer nov08
     ;;these next 2 lines are observatory specific
      if keyword_set(ctio4k) then dopenv=ctio4k_init(obsnm,'iod', 0.0, tag=tag, $
      	psfmod=psfmod, iss_obnm=obsnm, date=date, shft_style=shft_style) 
      if keyword_set(ctio) then dopenv=ctio_init(obsnm,'iod', 0.0, iss_obnm=obsnm) 
      if keyword_set(keck) then dopenv=keck_init(obsnm,'iod', 0.0, iss_obnm=obsnm, tag=tag) 
      if keyword_set(lick) then dopenv=lick_init(obsnm,'iod', 0.0, iss_obnm=obsnm, tag=tag) 
      restore,dopenv.bary_file  ;bcat
   	  x = where(strt(bcat.obsnm) eq 'r'+strmid(obsnm, 1, strlen(obsnm)-1), nx) 
      if nx le 0 then begin
	      x = where(strt(bcat.obsnm) eq strmid(obsnm, 1, strlen(obsnm)-1), nx) 
      endif

	if ~keyword_set(avg) then avg=0
      if nx gt 0 then begin
         bcat=bcat[x[0]]
         objnm = bcat.objnm
         found=findfile(dopenv.files_dir+'cd'+tag+'b*'+obsnm, count=count)
         if count eq 1 then print, found
         if count eq 0 then begin
		    print,obsnm
            dop_main, dopenv.obj_nm, dopenv, obsnm=obsnm, verbose=verbose, yrmo=yrmo, avg=avg, $
                              pass=1, tag=tag, /iod_soln, demo=demo,  $
                              vdavgnm=vdavgnm, cdnear_name=cdnear_name 
            dop_main, dopenv.obj_nm, dopenv, obsnm=obsnm, pass=2, tag=tag, /iod_soln 
            print, 'Done: '+found
      		if keyword_set(ctio4k) then dop_create_ipcf,tag=tag, observatory='ctio4k'
         endif
      endif
      if nx eq 0 then print,'Iodine obs not found to run'
end

