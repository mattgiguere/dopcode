pro  rayclean,obsnm,obstack,star,nocoadd=nocoadd,diff=diff,auto=auto, date=date

;obsnm [input]  string array    observation names used to be "ray cleaned
;                               and coadded to form template star
;star  [output] floating array  coadded stellar template
;        
;Created PB, Nov 27, 1995

;ordr = indgen(26)+28    ;operate on these orders
path='/mir1/fitspec_tst/'
ordr = indgen(16)+38     ;operate on LICK STANDARD orders
ordrstnd = 45            ;standard order
;mdtck = 2.       ;use this if coadding 5 or more spectra
mdtck = 3.5       ;use this if coadding 3 spectra
;mdtck = 1000
;mdtck = 2000
if strmid(obsnm(0),0,2) eq 'rk' then begin ;PB KLUDGE for Keck HR 493
;   ordr = indgen(10)  ;for pre-rk5
;   ordrstnd = 5       ;for pre-rk5 standard order for Keck rk2
   ordr = indgen(10)+21
   ordrstnd = 25                           ;standard order for Keck rk2
   mdtck=3000
   path='/mir3/fitspec/'
endif   ;Keck observation
if strmid(obsnm(0),0,2) eq 'rj' then begin ;
   ordr = indgen(15)
   ordrstnd = 3                           ;standard order for Keck rk2
   mdtck=3000
   path='/tous/mir3/fitspec/'
endif   ;Keck observation
if strmid(obsnm(0),0,2) eq 'rq' then begin ;CTIO for rq 
   ordr = indgen(18)+12
   ordrstnd = 22                           ;standard order 
   mdtck=1000
   path='/mir7/fitspec/'
endif   ;CTIO observation
if strmid(obsnm(0),0,4) eq 'achi' then begin ;CTIO for rq 
   ordr = indgen(18)+13
   ordrstnd = 22                           ;standard order 
   mdtck=1000
   path='/tous/mir7/fitspec/'+date+'/'
endif   ;CTIO observation
if strmid(obsnm(0),0,2) eq 'rc' then begin ;CTIO4k  
   ordr = indgen(18)+12
   ordrstnd = 22                           ;standard order 
   mdtck=1000
   path='/mir7/fitspec/'
endif   ;CTIO observation

nobs=n_elements(obsnm)
if nobs lt 2 then stop,'Less than 2 observations, unable to clean cosmic rays'
if nobs lt 3 and keyword_set(auto) then $
     stop,'Less than 3 observations, unable to automatically clean cosmic rays'
if nobs ge 3 then mdtck=3                  ;PB kludge, 22 Sep 98

ob=mrdfits(path+obsnm[0]+'.fits')
if n_elements(size(ob)) gt 5 then ob=reform(ob[1,*,*])
filter=ob*0+1                   ;perfect chip?
;rdsi,ob,obsnm(0)
mdcts=[median(ob(*,ordrstnd))]
PRINT,'median counts at order '+strtrim(ordrstnd,2)+' for '+obsnm(0)+'    '+strtrim(median(ob(*,ordrstnd)),2)
obstack = fltarr(n_elements(ob(*,0)),n_elements(ob(0,*)),nobs)
obstack(*,*,0) = ob
for n=1,nobs-1 do begin
   ob=mrdfits(path+obsnm[n]+'.fits')
;	rdsk,ob,path+obsnm[n],1
    if n_elements(size(ob)) gt 5 then ob=reform(ob[1,*,*])
   pfilt=ob*0+1                ;perfect chip?
   mdcts=[mdcts,median(ob(*,ordrstnd))]
   print,'Median counts at order '+strtrim(ordrstnd,2)+' for '+obsnm(n)+'    '+strtrim(median(ob(*,ordrstnd)),2)
   obstack(*,*,n)=ob
endfor

if nobs eq 2 then begin                         ; 2 observations, use RAYZAP
   print,'Cleaning with RAYZAP'
   ob=reform(obstack(*,*,0))
   dum=reform(obstack(*,*,1))
;   if n_elements(diff) ne 1 then diff=min([1400,median(ob(*,45))/10.])
   if n_elements(diff) ne 1 then diff=3.5
;   if diff lt 200 then diff=min([1400,median(ob(*,ordrstnd))/10.])
;   rayzap,ob,dum,pfilt,ordr,diff,/median
;   rayzap,ob,dum,pfilt,ordr,12000   ;for Keck template rk11.193, rk11.199
;   rayzap,ob,dum,pfilt,ordr,5000   ;for Keck template rk12 gl699
;   rayzap,ob,dum,pfilt,ordr,2000    ;for Keck template rk20 hd93745
   mincts=min([2000,0.25*median(mdcts)])  ;2000 DN or 1/4 of median exposure
   if n_elements(diff) eq 1 then if diff gt 50 then mincts=diff     ;PB 10Aug99
   rayzap,ob,dum,pfilt,ordr,mincts  ;for Keck template rk20 hd93745
   obstack(*,*,0) = ob
   obstack(*,*,1) = dum
endif else raystack,obstack,pfilt,ordr,mdtck=mdtck,auto=auto ; >2 observations, use RAYSTACK
;endif else raystack,obstack,pfilt,ordr,auto=auto  ; >2 observations, use RAYSTACK

;coadding the "cleaned" observations
if not keyword_set(nocoadd) then begin
   print,'Now Co-adding the observations!'
   dum=reform(obstack(*,*,0))
   for n=1,nobs-1 do begin
      coadd,dum,reform(obstack(*,*,n)),pfilt,ordr,sh,star  &  dum=star
   endfor
endif

return
end
