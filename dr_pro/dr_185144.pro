pro dr_185144, run=run, obnm=obnm, pass=pass, debug=debug, demo=demo, $
               jdmin=jdmin, jdmax=jdmax

  star='185144'
  if keyword_set(demo) then demo=1 else demo=0

  restore,'/mir3/bary/kbcvel.dat'
  xmatch=where(strmid(bcat.obsnm,0,2) eq 'rj' and bcat.objnm eq star $
               and bcat.jd gt jdmin and bcat.jd lt jdmax,num) 
print,num
  if num eq 0 then stop
  dopcount=0

  obs_id=bcat[xmatch].obsnm
print,obs_id
  for j=0,num-1 do begin 
     iss_nm='dsst185144ab_rj55.dat'
     iss_bc= 4082.8
     iss_obnm='rj55.1842'
     dopenv=keck_init(obs_id[j],iss_nm, iss_bc, iss_obnm=iss_obnm)
     tag='k'
     restore,dopenv.bary_file   ;bcat
     x=where(bcat.obsnm eq obs_id[j],nx)
     if nx gt 0 then begin
        bcat=bcat[x[0]]
        objnm = bcat.objnm
        found=findfile('/mir3/files/vd'+tag+objnm+'*'+obs_id[j], count=count)
        if count eq 1 then print, found
        ffits=findfile('/mir3/fitsfiles/'+obs_id[j]+'.fits',count=fcount)
        if count eq 0 and fcount eq 1 then begin  ;not run, but fitsfile exists
           if keyword_set(debug) then $
            dop_main, '185144',dopenv, obsnm=obs_id[j], pass=1, /demo, /verbose else $
               dop_main, '185144', dopenv, obsnm=obs_id[j], pass=1, tag=tag, demo=demo
           dop_main, '185144', dopenv, obsnm=obs_id[j], pass=2, tag=tag
           dopcount=dopcount+1    ;only count new files
        endif
     endif
     if dopcount gt 50 then stop
  endfor


end
