pro dr_obsa, obsnm=obsnm, tag=tag, pass=pass, debug=debug, demo=demo, date=date,$
	iodarr=iodarr

; EXAMPLE:
; dr_obsa, obsnm='aqa33.1059', tag='ad', /demo

; SEE CTIO4K_INIT for a description of the tags.
  star='128620'
  object=star
; if keyword_set(demo) then demo=1 else demo=0
  demo=1

     if tag eq 'j' then iss_nm='dsst128620b_rqa10.dat'  ; from JJ's morph code
     if tag eq 'b' then begin
     	iss_nm='dsst128620n_rqa10.dat'   ;dsst created Jan 24, 2010
     	iss_bc=-20540
     	iss_obnm='rqa10.6622'
     endif
;     if strmid(tag,0,1) eq 'a' or strmid(tag,0,1) eq 'b' then begin
;      ;80 chunk size
;     	iss_nm='dsst128620kd_rqa31.dat'
;     	iss_bc=21105.990
;     	iss_obnm='rqa31.127'
;     endif
     if strmid(tag,0,1) eq 'd' then begin
      ;240 chunk size
     	iss_nm='dsst128620dd_rqa31.dat'
     	iss_bc=21105.990
     	iss_obnm='rqa31.127'
     endif
     if strmid(tag,0,1) eq 'k' then begin
      ;80 chunk size
     	iss_nm='dsst128620kd_rqa31.dat'
     	iss_bc=21105.990
     	iss_obnm='rqa31.127'
     endif
     if strmid(tag,0,1) eq 'a' then begin  ; tag=ax2, adf, aff, afp
      ;80 chunk size
        restore,'/tous/mir7/bary/qbcvel.dat'  ;bcat
        x=where(bcat.obsnm eq 'chi120606.1125',nx)
     	iss_nm='dsst128620afm_achi120606.1125.dat'
     	iss_bc=bcat[x].bc
     	iss_obnm='achi120606.1125'
     	tmpl_dir='/tous/mir7/fitspec/120606/'  ; to find fits file for is_spec in dop_chunk_setup
        dpad=12
     endif
     if ~keyword_set(dpad) then dpad=12
     ;stop
     dopenv=ctio4k_init(obsnm,iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag, dpad=dpad, date=date)
     restore,dopenv.bary_file   ;bcat
;     x = where(bcat.obsnm eq 'r'+strt(strmid(obsnm, 1, strlen(obsnm)-1)), nx) 
     x = where(bcat.obsnm eq strt(strmid(obsnm, 1, strlen(obsnm)-1)), nx) 
     ;x=where(bcat.obsnm eq obsnm,nx)
;   help,/st,dopenv
 ;    stop
     if nx eq 0 then stop,'observation is not in qbcvel.dat' 
        found=findfile('/tous/mir7/files/vd'+tag+star+'*'+obsnm, count=count)
        if count eq 1 then print, found
        if count eq 0 then begin
           dop_main, '128620', dopenv, obsnm=obsnm, pass=1, tag=tag, demo=demo, $
           		verbose=verbose, tmpl_dir=tmpl_dir                   ;iodarr=iodarr, 
           dop_main, '128620', dopenv, obsnm=obsnm, pass=2, tag=tag
        endif

end
