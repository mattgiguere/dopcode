pro dr_obsb_loop, tag=tag, run=run, obmin=obmin, obmax=obmax, demo=demo, lm=lm, verbose=verbose

  star='128621'
  if keyword_set(demo) then demo=1 else demo=0
  ; default Levenberg-Marquardt fitting is with mpfit package
  if ~keyword_set(lm) then lm='mpf'
;  if ~keyword_set(obmin) then obmin=1000
;  if ~keyword_set(obmax) then obmax=9999

  restore,'/mir7/bary/qbcvel.dat'
  xmatch=where(strmid(bcat.obsnm,0,5) eq run and bcat.objnm eq '128621' and bcat.obtype eq 'o' $
               and fix(strmid(bcat.obsnm,6,4)) ge obmin and fix(strmid(bcat.obsnm,6,4)) le obmax,num)
  if num eq 0 then stop,'No matching spectra to analyze'
  dopcount=0
  obs_id=bcat[xmatch].obsnm

  for j=0,num-1 do begin 
     if tag eq 'b' then iss_nm='dsst128621n_rqa10.dat'   ;dsst created Apr 8, 2010 ;Jan 24, 2010
     iss_bc=-20576
     iss_obnm='rqa10.6662'
     dopenv=ctio_init(obs_id[j],iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag)
     restore,dopenv.bary_file   ;restore bary structure, "bcat"
     x=where(bcat.obsnm eq obs_id[j],nx)
     if nx gt 1 then stop,'duplicate entries - fix bary file'
     if nx eq 1 then begin
        bcat=bcat[x[0]]
        objnm = bcat.objnm
        found=findfile('/mir7/files/vd'+tag+objnm+'*'+obs_id[j], count=count)
        if count eq 1 then print, found
        ffits=findfile('/mir7/fitspec/'+obs_id[j]+'.fits',count=fcount)
        if count eq 0 and fcount eq 1 then begin  ;not run, but fitsfile exists
           dop_main, '128621',dopenv, obsnm=obs_id[j], pass=1, tag=tag, lm=lm, $
                     demo=demo, verbose=verbose
           dop_main, '128621', dopenv, obsnm=obs_id[j], pass=2, tag=tag, lm=lm
           dopcount=dopcount+1    ;only count new files
        endif
        if count eq 0 and fcount eq 0 then stop,'need to run obs, but fitsfile not found'
     endif
     print,'DOPCOUNT: ',dopcount,' OUT OF: ',num
  endfor


end
