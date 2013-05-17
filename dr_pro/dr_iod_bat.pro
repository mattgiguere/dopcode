pro dr_iod_bat, obsnm, tag=tag, demo=demo, verbose=verbose, ctio=ctio,keck=keck 

; fischer nov08
      iss_nm='iod'
      iss_bc=0.0
      iss_obnm='iod'

     ;;these next 2 lines are observatory specific
      if keyword_set(ctio) then dopenv=ctio_init(obsnm, iss_nm, iss_bc, iss_obnm=iss_obnm) 
      if keyword_set(keck) then dopenv=keck_init(obsnm, iss_nm, iss_bc, iss_obnm=iss_obnm) 
      if keyword_set(het) then dopenv=het_init(obsnm, iss_nm, iss_bc, iss_obnm=iss_obnm) 
      if ~keyword_set(tag) then tag='c'
      dop_main, dopenv.obj_nm, dopenv, obsnm=obsnm, $
                pass=1, tag=tag, /iod_soln, demo=demo 
      dop_main, dopenv.obj_nm, dopenv, obsnm=obsnm, pass=2, tag=tag, /iod_soln 
      if nx eq 0 then print,'Iodine obs not found to run'

end
