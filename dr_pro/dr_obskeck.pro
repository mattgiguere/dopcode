pro dr_obskeck, obnm=obnm, star=star, pass=pass, demo=demo

tag='js'  ; tag for ipcf and observations
dssttag='js'

restore,'/tous/mir3/bary/kbcvel.dat'

ff=file_search('/tous/mir3/files/dsst'+star+dssttag+'_rj162*',count=count)
if count eq 0 then stop,'no dsst found'
iss_nm = strmid(ff,17,strlen(ff)-17)
print,iss_nm
xx=where(bcat.objnm eq star and bcat.obtype eq 't' and strmid(bcat.obsnm,0,5) eq 'rj162',nxx)

; get info about the dsst
;     for j=0,nxx-1 do begin
;        print,j,bcat[xx[j]].obsnm,bcat[xx[j]].jd,format='(i4,a12,f15.5)'
;     endfor
;     ans=0
;     read,'which templ obs should be used? (e.g., 2): ',ans
;     xx=xx[ans]
  iss_obnm = bcat[xx[2]].obsnm
  x=where(bcat.obsnm eq iss_obnm,nx)  ; to get bc of dsst
  if nx eq 0 then stop,'no dsst observation found' else iss_bc=bcat[x].bc
; done with dsst info

  if keyword_set(demo) then demo=1 else demo=0
  ; default Levenberg-Marquardt fitting is with mpfit package
  xmatch=where(bcat.obsnm eq obnm,num)
  if num eq 0 then stop,'No stellar spectra to analyze'
  bcat=bcat[xmatch]
  dopenv=keck_init(obnm, iss_nm, iss_bc, iss_obnm=iss_obnm, tag=tag)
  restore,dopenv.bary_file      ;restore bary structure, "bcat"

  found=findfile('/tous/mir3/files/vd'+tag+'*'+obnm, count=count)
  if count eq 1 then print, found
  ffits=findfile('/tous/mir3/fitspec/'+obnm+'.fits',count=fcount)
  if count eq 0 and fcount eq 1 then begin ;not run, but fitsfile exists
     dop_main, star ,dopenv, obsnm=obnm, pass=1, tag=tag, lm=mpf, $
               demo=1, verbose=verbose
     dop_main, star, dopenv, obsnm=obnm, pass=2, tag=tag, lm=mpf, demo=0
  endif
  if count eq 0 and fcount eq 0 then stop,'need to run obs, but fitsfile not found'

end

