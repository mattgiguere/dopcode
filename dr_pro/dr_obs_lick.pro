pro dr_obs_lick, star, obsnm=obsnm, pass=pass, debug=debug, demo=demo

; fischer 17nov2011,  9feb2013
; dr_obs_lick_loop, '161797', obsnm='ri69.86', iss_nm='dsst161797fib_ri68.dat', $
;                    tag='dff', vdavg='vddff_iod.ri69.n1.avg'

  if keyword_set(demo) then demo=1 else demo=0
	tag='n'
	iss_nm='dssthip96459n_rs03.52.dat'
	iss_bc= -2430.254
	iss_obnm='rs03.52'
	
     dopenv=lick_init(obsnm,iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag)
     restore,dopenv.bary_file   ;bcat
     x=where(bcat.obsnm eq obsnm,nx)
     if nx eq 0 then stop,'observation is not in qbcvel.dat' 
        found=findfile('/mir1/files/vd'+tag+star+'*'+obsnm, count=count)
        if count eq 1 then print, found
        if count eq 0 then begin
           dop_main, star, dopenv, obsnm=obsnm, pass=1, tag=tag, demo=demo, verbose=verbose
           dop_main, star, dopenv, obsnm=obsnm, pass=2, tag=tag
        endif
   

end
