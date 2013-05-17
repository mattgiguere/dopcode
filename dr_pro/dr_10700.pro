pro dr_10700, obsnm=obsnm,run=run, tag=tag, pass=pass, debug=debug, demo=demo, date=date,$
	iodarr=iodarr, solarmorph=solarmorph, jansson=jansson, verbose=verbose

; EXAMPLE:
; dr_10700,obsnm='achi120710.1148',tag='am',date='120710',/demo
; dr_10700,obsnm='rchi111104.1255',tag='am',date='111104',/demo
; dr_10700,run='120710',tag='adg',date='120710'

; SEE CTIO4K_INIT for a description of the tags.
  star='10700'
  object=star
 if keyword_set(demo) then demo=1 else demo=0

     if strmid(tag,0,1) eq 'a' then begin  ; tag=ax2, adf, aff, afp, adg
        restore,'/tous/mir7/bary/qbcvel.dat'  ;bcat
        x=where(bcat.obsnm eq 'chi120710.1141',nx)
        iss_nm = 'dsst10700adg_achi120710.1141.dat'
     	iss_bc=bcat[x].bc
print,iss_bc
     	iss_obnm='achi120710.1141'
     	tmpl_dir='/tous/mir7/fitspec/120710/'  ; to find fits file for
        dpad=12
     endif

     if strmid(tag,0,1) eq 'k' then begin  ; tag=ax2, adf, aff, afp, am
      ;80 chunk size
        restore,'/tous/mir7/bary/qbcvel.dat'  ;bcat
        x=where(bcat.obsnm eq 'rchi111102.1272',nx)
        iss_nm = 'dsst10700kd1jan_rchi111102.1272.dat'
     	iss_bc=bcat[x].bc
print,iss_bc
     	iss_obnm='rchi111102.1272'
     	tmpl_dir='/tous/mir7/fitspec/'  ; to find fits file for is_spec in dop_chunk_setup
        dpad=12
     endif
     
     if ~keyword_set(dpad) then dpad=12
     ;stop
     dopenv=ctio4k_init(obsnm,iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag, dpad=dpad, date=date)
     restore,dopenv.bary_file   ;bcat
;     x = where(bcat.obsnm eq 'r'+strt(strmid(obsnm, 1, strlen(obsnm)-1)), nx) 
     if strmid(obsnm,0,1) eq 'a' then x = where(bcat.obsnm eq strmid(obsnm,1,strlen(obsnm)-1), nx) $
     else x = where(bcat.obsnm eq obsnm, nx) 

     if nx eq 0 then stop,'observation is not in qbcvel.dat' 
        found=findfile('/tous/mir7/files_df/vd'+tag+star+'*'+obsnm, count=count)
        if count eq 1 then print, found
        if count eq 0 then begin
           dop_main, '10700', dopenv, obsnm=obsnm, pass=1, tag=tag, demo=demo, $
           		verbose=verbose, tmpl_dir=tmpl_dir, vdtag=tag                   ;iodarr=iodarr, 
           dop_main, '10700', dopenv, obsnm=obsnm, pass=2, tag=tag, vdtag=tag
        endif

end
