pro dr_star, star=star, obsnm=obsnm, tag=tag, date=date, yrmo=yrmo, cdnear_name=cdnear_name, $
			 dsst_nm=dsst_nm, vdtag=vdtag, avg=avg, vdavgnm=vdavgnm, demo=demo, $
			 psfmod=psfmod, shft_style=shft_style, verbose=verbose


;Fischer 2013 Jan
; called by dr_obs to run all observations in a night 
; or run alone to analyze a single observation

; dr_star, star='10700', obsnm='achi120710.1148', tag='adg',date='120710', $
;          yrmo=yrmo, dsst_nm='dsst10700adg_achi120710.1141.dat'

; SEE CTIO4K_INIT for a description of the tags.
 if ~keyword_set(demo) then demo=0
 ff=file_search('/tous/mir7/files_df/'+dsst_nm,count=count)
 ;print,count, '  ', ff
 if count eq 0 then stop, 'no dsst found'
 
 d=strpos(dsst_nm,'chi')
 nm=strmid(dsst_nm,d,strlen(dsst_nm)-d-4) 
 d2=strpos(nm,'.')
 dsst_date = strmid(nm,3,6)

        restore,'/tous/mir7/bary/qbcvel.dat'  ;bcat
        x=where(bcat.obsnm eq nm,nx)
        iss_nm = dsst_nm
     	iss_bc=bcat[x].bc
     	iss_obnm='a'+nm
     	tmpl_dir='/tous/mir7/fitspec/'+dsst_date+'/'  ; to find fits file for dsst
        dpad=12
     
     dopenv=ctio4k_init(obsnm,iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag, dpad=dpad, $
     	psfmod=psfmod, date=date, avg=avg, shft_style=shft_style)
     restore,dopenv.bary_file   ;bcat
     if strmid(obsnm,0,1) eq 'a' then x = where(bcat.obsnm eq strmid(obsnm,1,strlen(obsnm)-1), nx) $
     else x = where(bcat.obsnm eq obsnm, nx) 
     if nx eq 0 then x=where(bcat.obsnm eq 'r'+strmid(obsnm,1,strlen(obsnm)-1),nx)
     if nx eq 0 then stop,'observation is not in qbcvel.dat' 
	 if ~keyword_set(vdtag) then vdtag=tag

     found=findfile('/tous/mir7/files_df/vd'+vdtag+star+'*'+obsnm, count=count)
     if count eq 1 then print, found
     if count eq 0 then begin
        dop_main, star, dopenv, obsnm=obsnm, pass=1, yrmo=yrmo, cdnear_name=cdnear_name, $
        	tag=tag, demo=demo, verbose=verbose, tmpl_dir=tmpl_dir, vdtag=vdtag, avg=avg, $
        	vdavgnm=vdavgnm
        dop_main, star, dopenv, obsnm=obsnm, pass=2, yrmo=yrmo, tag=tag, vdtag=vdtag
     endif

end
