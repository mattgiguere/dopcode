pro dr_obsa_bat, obnm=obnm, pass=pass, debug=debug, demo=demo

  if keyword_set(demo) then demo=1 else demo=0

     iss_nm='dsst128620n_rqa10.dat'
     iss_bc=-20540
     iss_obnm='rqa10.6622'
     dopenv=ctio_init(obnm,iss_nm, iss_bc, iss_obnm=iss_obnm)
     tag='c'
     restore,dopenv.bary_file   ;bcat
     x=where(bcat.obsnm eq obnm,nx)
     if nx gt 0 then begin
        bcat=bcat[x[0]]
        objnm = bcat.objnm
        dop_main, '128620', dopenv, obsnm=obnm, pass=1, tag=tag, demo=demo
        dop_main, '128620', dopenv, obsnm=obnm, pass=2, tag=tag
     endif

end
