PRO DOP_RAYCLEAN,  obs_arr, dopenv=dopenv, observatory=observatory,obstack=obstack,star=star,nocoadd=nocoadd,diff=diff,auto=auto

; PROCEDURES CALLED: 
; rayzap
; coadd -> chicorl, shfour, idlcorl
; raystack
; rayslave
;
; INPUT
; OBS_ARR: array of N observation names to be raycleaned and co-added
; DOPENV: structure containing info about the observations
; OBSERVATORY: ctio4k, lick, keck
;
;OUTPUT: 
; OBSTACK: the array of N spectra
; STAR: variable name for the co-added float array
;        
; OPTIONS
; NOCOADD: 
; DIFF: 
; AUTO:  
; 
;Created PB, Nov 27, 1995
;cleaned up and modified Fischer Jan 30, 2012 

path=dopenv.obs_dir+'/'   ; fits format
;if n_elements(obs_arr ge 5) then mdtck = 2. 
;if n_elements(obs_arr lt 5) then 

if observatory eq 'lick' then begin
   ordr = indgen(16)+38         ;operate on LICK STANDARD orders
   ordrstnd = 45                ;standard order
   mdtck = 3.5
endif
if observatory eq 'keck' then begin 
;   ordr = indgen(10)  ;for pre-rk5
;   ordrstnd = 5       ;for pre-rk5 standard order for Keck rk2
   ordr = indgen(10)+21
   ordrstnd = 25                           ;standard order for Keck rk2
   mdtck=3000
endif   ;Keck observation
if observatory eq 'ctio4k' then begin 
   ordr = indgen(18)+12
   ordrstnd = 22                           ;standard order 
   mdtck=1000
endif   ;CTIO observation

nobs=n_elements(obs_arr)
if nobs lt 2 then stop,'Less than 2 observations, unable to clean cosmic rays'
if nobs lt 3 and keyword_set(auto) then $
     stop,'Less than 3 observations, unable to automatically clean cosmic rays'
if nobs ge 3 then mdtck=3                  ;PB kludge, 22 Sep 98

ob=readfits(path+obs_arr[0]+'.fits')
if n_elements(size(ob)) gt 5 then ob=reform(ob[1,*,*])
filter=ob*0+1                   ;perfect chip?

mdcts=[median(ob[*,ordrstnd])]
PRINT,'median counts at order '+strtrim(ordrstnd,2)+' for '+obs_arr[0]+'    '+strtrim(median(ob[*,ordrstnd]),2)
obstack = fltarr(n_elements(ob[*,0]),n_elements(ob[0,*]),nobs)
obstack[*,*,0] = ob
for n=1,nobs-1 do begin
   ob=readfits(path+obs_arr[n]+'.fits')
    if n_elements(size(ob)) gt 5 then ob=reform(ob[1,*,*])
   pfilt=ob*0+1                ;perfect chip?
   mdcts=[mdcts,median(ob[*,ordrstnd])]
   print,'Median counts at order '+strtrim(ordrstnd,2)+' for '+obs_arr[n]+'    '+strtrim(median(ob[*,ordrstnd]),2)
   obstack[*,*,n]=ob
endfor

if nobs eq 2 then begin                         ; 2 observations, use RAYZAP
   print,'Cleaning with RAYZAP'
   ob=reform(obstack[*,*,0])
   dum=reform(obstack[*,*,1])
   mincts=min([2000,0.25*median(mdcts)])    ;2000 DN or 1/4 of median exposure
   rayzap,ob,dum,pfilt,ordr,mincts
   obstack[*,*,0] = ob
   obstack[*,*,1] = dum
endif
if nobs gt 2 then raystack,obstack,pfilt,ordr,mdtck=mdtck,auto=auto 

;coadding the "cleaned" observations
if ~keyword_set(nocoadd) then begin
   print,'Co-adding the observations'
   dum=reform(obstack[*,*,0])
   for n=1,nobs-1 do begin
      coadd,dum,reform(obstack[*,*,n]),pfilt,ordr,sh,star  
      dum=star
   endfor
endif

return
end
