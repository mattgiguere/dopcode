pro dr_obskeck_loop, star=star,iss_nm=iss_nm, iss_obnm=iss_obnm,$
                 tag=tag, dsst_tag=dsst_tag, jdmin=jdmin, demo=demo, $
                 lm=lm, verbose=verbose,crosscorl=crosscorl 

; dr_obs_loop,star='HIP75458',iss_nm='dssthip75458df_rs17.dat',$
;  iss_obnm='rs17.124', tag='b',jdmin=11000.
  restore,'/mir3/bary/kbcvel.dat'
  xx=where(bcat.objnm eq star and bcat.obtype eq 't',nxx)
  if nxx eq 0 then stop, 'No matches in bcvel'
  if nxx gt 1 then begin 
     for j=0,nxx-1 do begin
        print,j,bcat[xx[j]].obsnm,bcat[xx[j]].jd,format='(i4,a12,f15.5)'
     endfor
     ans=0
     read,'which templ obs should be used? (e.g., 2): ',ans
     xx=xx[ans]
  endif
  iss_obnm = bcat[xx].obsnm
  if ~keyword_set(dsst_tag) or ~keyword_set(iss_nm) then begin
     ans='' & read,'the default dsst_tag is fn, is this OK? (y/n) ',ans
     if ans eq 'y' then dsst_tag = 'fn'
     if ans eq 'n' then begin
        dsst_tag=''
        read,'enter the correct dsst_tag: ',dsst_tag
     endif
  endif
  iss_nm = 'dsst'+star+dsst_tag+'_'+strmid(bcat[xx].obsnm,0,4)+'.dat'
  iffnm=file_search('/mir3/files/'+iss_nm,count=count)   
  if count eq 0 then $
  iss_nm = 'dsst'+star+dsst_tag+'_'+strmid(bcat[xx].obsnm,0,5)+'.dat' 
  iffnm=file_search('/mir3/files/'+iss_nm,count=count)  
  if count eq 0 then stop, 'iss not found'

  if ~keyword_set(jdmin) then jdmin=15000.

  x=where(bcat.obsnm eq iss_obnm,nx)  ; to get bc of dsst
  if nx eq 0 then stop,'no dsst found' else iss_bc=bcat[x].bc
  
  if keyword_set(demo) then demo=1 else demo=0
  ; default Levenberg-Marquardt fitting is with mpfit package
  if ~keyword_set(lm) then lm='mpf'
  if ~keyword_set(crosscorl) then crosscorl=0 

  xmatch=where(bcat.jd gt jdmin and bcat.objnm eq star and bcat.obtype eq 'o',num)

  if num eq 0 then stop,'No stellar spectra to analyze'
  dopcount=0
  obs_id=bcat[xmatch].obsnm

  for j=0,num-1 do begin 
;stop
     dopenv=keck_init(obs_id[j], iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag)
     restore,dopenv.bary_file   ;restore bary structure, "bcat"
     x=where(bcat.obsnm eq obs_id[j],nx)
     if nx gt 1 then stop,'duplicate entries - fix bary file'
     if nx eq 1 then begin
        bcat=bcat[x[0]]
        objnm = bcat.objnm
        found=file_search('/mir3/files/vd'+tag+objnm+'*'+obs_id[j], count=count)
        if count eq 1 then print, found
        ffits=file_search('/mir3/fitspec/'+obs_id[j]+'.fits',count=fcount)
        if count eq 0 and fcount eq 1 then begin  ;not run, but fitsfile exists
           dop_main, star ,dopenv, obsnm=obs_id[j], pass=1, tag=tag, lm=lm, $
                     demo=demo, verbose=verbose, crosscorl=crosscorl
           dop_main, star, dopenv, obsnm=obs_id[j], pass=2, tag=tag, lm=lm
           dopcount=dopcount+1    ;only count new files
        endif
        if count eq 0 and fcount eq 0 then stop,'need to run obs, but fitsfile not found'
     endif
     print,'DOPCOUNT: ',dopcount,' OUT OF: ',num
  endfor


end
