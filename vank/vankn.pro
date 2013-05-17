pro vankn,strnm,tmptp,lbl,cf1,cf3,mct=mct,lazer=lazer,maxchi=maxchi,$
noclean=noclean,emu=emu,cf=cfout,ps=ps,title=title,$
outfile=outfile,vdarr=vdarr, reject=reject,dircf=dircf,fit_key=fit_key,$
doptest=doptest,nopatch=nopatch

;fit_key = 1 for stars with lick templates
;fit_key = 2 for stars with keck templates
;fit_key = 3 for new reduction

if ~keyword_set(tmptp) then tmptp='f'
if ~keyword_set(lbl) then lbl='vdf'
if ~keyword_set(fit_key) then fit_key=1
if ~keyword_set(doptest) then doptest=1

;vank,'75732','daf','daf',cf1,cf3 

; cf1  output structure   all velocities
; cf2  output structure    nightly corrected version of cf1
; cf3  output structure    "Raw Counts" corrected version of cf1

; mct  keyword integer     minimum acceptable counts,
;			     overrides the hardwired rules

;INPUTS
;strnm  (string)  '458'    star_name
;tmptp  (string)  'rk12'   template_tape
;lbl    (string)  'vdv'    VD label

;;When using /emu and test other than 'em1', do this e.g.:
;; vank,'emu','ste',/emu
;; When the VD files were named vdlm_em0.19, I used
;; vank,'lm','em0','vd',/emu


mk_cf,strnm

;vstdsk = '/further/paul/idle/kdop/kvel/vst'
vstdsk = '~fischer/vank/'
cfdsk = '/mir1/dop/cf'
fildsk = '/mir1/files/'+lbl+strnm+'_'
;fildsk = '/mir1/files/vda_etoile/'+lbl+strnm+'_'
kecklist = '/mir3/files/keck_st.dat'
;rawfile  = '/mir3/logsheets/rawcts.dat'

;restore,kecklist                ;get Keck List Data Structure
;keckst=dum                      ;rename Keck List Data Structure 

;HD Star Catalog
catalog = 'HD '  
;cfnm    = cfdsk  + strupcase(strtrim(strnm,2))+'_'+tmptp+'.dat'
if keyword_set(emu) then begin
    cfnm = cfdsk + 'EMU_' + tmptp + '.dat'
    maxchi = 40
endif else cfnm    = cfdsk  + strlowcase(strmid(strnm,0,3)) + $
                     strtrim(strupcase(strmid(strnm,3,20)),2) +'_'+tmptp+'.dat'
if keyword_set(dircf) then cfnm = dircf

if strmid(strnm,0,2) eq '38' then $ ;special kludge for HD38A and HD38B
  cfnm = cfdsk  + strupcase(strtrim(strnm,2))+'_'+tmptp+'.dat'
vstnm   = vstdsk + 'vst'+strlowcase(strtrim(strnm,2))+'.dat'
vnm     = strnm
if 1-keyword_set(emu) then mncts   = 3000
;if strupcase(strtrim(strnm,2))  eq '219B' then mncts=1500     ;faintest HR star

if strupcase(strtrim(strnm,2))  eq 'GL876' then mncts=1500 ;faintest HR star

if strupcase(strmid(strnm,0,3)) eq 'HIP' then mncts=1000 
if strupcase(strnm) eq 'GL905' then mncts=1000 
if strupcase(strmid(strnm,0,2)) eq 'BD' then mncts=1000 
if strupcase(strmid(strnm,0,2)) eq 'G2' then mncts=1000 
if strupcase(strtrim(strnm,2))  eq '75732B' then mncts=1000 ;faint star

;Don't print HD for GLIESE or SAO or HIP or BD stars
if strupcase(strmid(strnm,0,1)) eq 'G' or $
  strupcase(strmid(strnm,0,1)) eq 'S' or $
  strupcase(strmid(strnm,0,3)) eq 'HIP' or $
  strupcase(strmid(strnm,0,2)) eq 'BD' then catalog=''

if n_elements(mct) eq 1 then if mct gt 0 then begin
    print,'Using input value for minimum acceptable photons in exposure:  '+strtrim(mct,2)
    mncts=mct
endif

dum=check_file(cfnm)

if dum then cmrestore,cfnm else begin
    cfnm = '/mir1/dop/cf'+strmid(strnm,0,3) + $
           strtrim(strlowcase(strmid(strnm,3,20)),2) +'_'+tmptp+'.dat'
    if check_file(cfnm) then restore,cfnm else begin
        message,'CF file ('+cfnm+' not located'
        return
    endelse 
endelse

cff=cf

if keyword_set(reject) then begin
    match, cff.obnm, reject, a, b
    if a[0] ge 0 then remove, a, cff
endif

;Toss the bad apples
cff=cff(where(cff.obnm ne 'rk12.400')) ;HD 105113 bad point chi-sq > 25!
cff=cff(where(cff.obnm ne 'rk11.63')) ;HD 139323 bad point chi-sq > 8!
cff=cff(where(cff.obnm ne 'rk11.151')) ;HD 172310 wrong star? SB2?, 12jul98
cff=cff(where(cff.obnm ne 'rk12.166')) ;HD 172310 wrong star? SB2?, 12jul98
cff=cff(where(cff.obnm ne 'rk12.353')) ;HR 7504 mididentified as HR 7503
cff=cff(where(cff.obnm ne 'rk12.520')) ;HR 7504 mididentified as HR 7503
cff=cff(where(cff.obnm ne 'rk11.103')) ;HD 63077 wrong star
cff=cff(where(cff.obnm ne 'rk13.109')) ;HD 218209 very noisy, 12jul98
cff=cff(where(cff.obnm ne 'rk13.112')) ;HD 1326B, fit=6.8
cff=cff(where(cff.obnm ne 'rk13.113')) ;HD 1326
cff=cff(where(cff.obnm ne 'rk26.941')) ;HD 210277, Cochran & Lissauer 900 m/s
cff=cff(where(cff.obnm ne 'rk27.20')) ;HD 217107 saturated
cff=cff(where(cff.obnm ne 'rk14.152')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk14.153')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk14.156')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk14.157')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk14.158')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk14.159')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk15.210')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk15.211')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk15.212')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk15.213')) ;HD 10476 test
cff=cff(where(cff.obnm ne 'rk11.104')) ;HD 67458 high cts, errvel 2x mean, 22mar01
cff=cff(where(cff.obnm ne 'rk11.106')) ;HD 71334 high cts, errvel 2x mean, 22mar01
cff=cff(where(cff.obnm ne 'rk12.493')) ;HD 95735 errvel twice mean, mdchi > 3.5 22mar01
cff=cff(where(cff.obnm ne 'rk5.90')) ;HD 135101A, mdchi > 4 22mar01
cff=cff(where(cff.obnm ne 'rk32.212')) ;HD 219834B, mdchi~3, errvel 2x mean 22mar01
cff=cff(where(cff.obnm ne 'rk29.292')) ;HD 363954, ??? 22mar01
cff=cff(where(cff.obnm ne 'rk58.273')) ;Lost power to I2 cell
cff=cff(where(cff.obnm ne 'rk58.274')) ;Lost power to I2 cell
cff=cff(where(cff.obnm ne 'rk58.275')) ;Lost power to I2 cell
cff=cff(where(cff.obnm ne 'rk58.276')) ;Lost power to I2 cell
cff=cff(where(cff.obnm ne 'rk58.277')) ;Lost power to I2 cell

dd=cff.jd
mindd=min(cff.jd)  
if mindd gt 10200 and mindd lt 2440000 then begin
    mindd=mindd-10200.
    dd=dd-10200.
endif
;if mindd gt 2440000 then begin
;   mindd=mindd-2440000.
;   dd=dd-2440000.
;endif

tape = strmid(cff[0].obnm,0,2)
case 1 of
    tape eq 'rk': begin
        ordr=[21,30]  
        pixr=[70,1950]
    end
    tape eq 'rf' or tape eq 'rg' or tape eq 'ri': begin
        ordr=[38,53]
        pixr=[50,1770]
    end
    tape eq 'rj' : begin
        ordr=[0, 14]
        pixr=[50, 3970]
    end
    tape eq 'em' and cff[0].jd ge 2.44d6+13480: begin
        ordr=[0, 14]
        pixr=[50, 3970]
    end
    tape eq 'em' and cff[0].jd lt 2.44d6+13480: begin
        ;ordr=[21,30]  
        ;pixr=[70,1950]

        ; Lick
        ordr=[38,58]  
        pixr=[50,1770]
    end
    tape eq 'rb' or tape eq 'rd' or tape eq 're': begin
        ordr=[38,53]
        pixr=[50,1770]
    end
    tape eq 'rs' or tape eq 'rt': begin
        ordr=[38,53]
        pixr=[50,1770]
    end
    tape eq 'ra' or tape eq 'rh' or tape eq 'rc': begin
	tapestring='?'
    end
    else: stop
endcase

if n_elements(cff) gt 1 then begin
    if n_elements(cff) eq 2 then begin
        restore,fildsk+cff(0).obnm  &   vd0=vd
        restore,fildsk+cff(1).obnm  &   vd1=vd
        cf1=cff
        dumfit=3.*median([vd0.fit1,vd1.fit1])
        dumwt=0.5*median(vd0.weight)
        ind=where(vd0.fit1 lt dumfit and vd1.fit1 lt dumfit and $
                  vd0.weight gt dumwt and $
                  vd0.pixt gt pixr(0) and vd0.pixt lt pixr(1) and $
                  vd0.ordt ge ordr(0) and vd0.ordt le ordr(1),nchunk)
        cf1(0).mnvel=mean(vd0(ind).vel)  &  cf1(0).mdvel=median(vd0(ind).vel)
        cf1(1).mnvel=mean(vd1(ind).vel)  &  cf1(1).mdvel=median(vd1(ind).vel)
        cf1(0).errvel=stdev(vd0(ind).vel)/sqrt(nchunk-1)
        cf1(1).errvel=stdev(vd1(ind).vel)/sqrt(nchunk-1)
        cf1(0).mdchi=median(vd0(ind).fit1)
        cf1(1).mdchi=median(vd1(ind).fit1)
        cf1(0).nchunk=nchunk
        cf1(1).nchunk=nchunk
        cf1(0).cts=median(vd0(ind).cts)
        cf1(1).cts=median(vd1(ind).cts)
        for n=0,19 do begin
            cf1(0).mnpar(n)=mean(vd0(ind).iparam(n))
            cf1(1).mnpar(n)=mean(vd1(ind).iparam(n))
            cf1(0).mdpar(n)=median(vd0(ind).iparam(n))
            cf1(1).mdpar(n)=median(vd1(ind).iparam(n))
        endfor
        cf1.jd=cf1.jd-2440000.
    endif else begin

if keyword_set(fit_key) then vel,vnm,cff,lbl,cf1,vdarr,ordr=ordr,pixr=pixr,$
  mincts=mncts,maxchi=maxchi,ifit=fit_key,noclean=noclean,doptest=doptest else $
	     vel,vnm,cff,lbl,cf1,vdarr,ordr=ordr,pixr=pixr,mincts=mncts,maxchi=maxchi,$
	     ifit=ifit,noclean=noclean,doptest=doptest
	end ;else
endif else begin
    print,'Less than two observations'
    return
endelse

;velplot,cf1,catalog+strupcase(strtrim(strnm,2)),0.6,d1,v1,e1,errcut=2.5,/yrs

;Get "Raw Counts" for correction
;  note: only "cf.sp2", the 99% level is used in correction
;             "cf.sp1" is used to store H&K Svalue
;restore,'rawcts.dat'
;restore,rawfile
;for n=0,n_elements(cf1)-1 do begin
;    qq=where(rawcts.obnm eq cf1(n).obnm,nqq)
;    if nqq eq 1 then cf1(n).sp1=rawcts(qq).cts98 ;initially store 98%
;    if nqq eq 1 then cf1(n).sp2=rawcts(qq).cts99 ;use only 99% for correction
;    if nqq ne 1 then cf1(n).sp1=0
;    if nqq ne 1 then cf1(n).sp2=0
;    if nqq ne 1 then print,'No Raw Counts for: '+cf1(n).obnm
;    if 1-keyword_set(emu) then begin
;;        hk,cf1(n).obnm,Sv       ;get the H&K Svalue
;;        cf1(n).sp1=Sv           ; store it in "cf.sp1"
;    endif
;endfor

;Nightly Correction, should be tried again, PB -- 8 Sept 99
cf5=cf1
;print,'Making **Primitive** Nightly Correction'
;knave,toss=vnm                             ;exclude post-fix
;rascii,knave,4,navdsk;   ,skip=1
;for n=0,n_elements(cf5)-1 do begin
;   nnd=minloc(abs(reform(knave(0,*)) - cf5(n).jd),/first) ;make correction
;   if abs(cf5(n).jd-knave(0,nnd)) lt 0.5 then cf5(n).mnvel=cf5(n).mnvel-knave(1,nnd) $
;   else print,'No nightly correction for the night of: '+strtrim(cf5(n).jd,2)
;endfor
;velplot,cf5,'CORRECTED: '+catalog+strupcase(strtrim(strnm,2)),0.6,d1,v1,e1,errcut=2.5,/yrs

;Apply Sine correction 
;print,'Making **Primitive** S/N Correction'
;sp=[10.723411,383.61427,3332.0877,-0.37066408]
;qq=where(cf1.jd lt 11413,nqq)    ;jd 11413 rk33, end of prefix
;qq=indgen(n_elements(cf1))      ;apply correction to all observations
;if nqq gt 0 then begin
;    x=sqrt(cf1(qq).cts)
;    yy=sincon(x,sp)
;    cf5(qq).mnvel=cf5(qq).mnvel-yy
;endif
;velplot,cf5,'SINE CORRECTED: '+catalog+strupcase(strtrim(strnm,2)),0.6,d1,v1,e1,errcut=2.5,/yrs

;Get B-V value
bv=0.                  

; dave - I don't have this. Should I?
;qq=first_el(where(keckst.name eq strtrim(strnm,2),nqq)) 
;if nqq eq 1 then if qq ge 0 then bv=keckst(qq).bv $
;                 else print,strnm+': Not Found in Keck List Data Structure' 

;Apply 11-th order polynomial correction to "Raw Counts"
cf3=cf1
;if bv lt 1.2 then begin   ;Don't apply corrections to late K and M stars 
;print,'Making **RAW Counts** S/N Correction'
;  11th order polynomial "correction" coefficents
;g11 = [-168.03319,1848.3498,-8358.2410,21065.190,-32578.292,32476.605]
;g11 = [g11,-21476.747,9513.7542,-2790.1381,519.48708,-55.567260,2.5984175]
;k11 = [37.249823,-518.41531,2634.3063,-6504.2037,9463.9504,-8956.8822]
;k11 = [k11,5741.9249,-2511.5797,736.70768,-138.27175,14.975206,-0.71060957]
;m11 = [291.32463,-3202.6396,14399.051,-34040.579,46120.661,-35073.525]
;m11 = [m11,10930.054,4325.1796,-5689.4488,2373.4906,-472.10427,37.509623]
;   qq=where(cf1.jd lt 11413,nqq)    ;jd 11413 rk33, end of prefix
;qq=indgen(n_elements(cf1))      ;apply correction to all observations
;if nqq gt 0 then begin
;    x=cf3(qq).sp2/1.e4          ; polynomial expects counts over 10,000
;    plcof=g11
;    if bv gt 0.8 then plcof=k11
;    if bv gt 1.2 then plcof=m11
;    yy=poly_fat(x,plcof)
;; Only apply correction for raw counts between 2,500 or 40,000 
;    i=where(x lt 2500/1.e4 or x gt 4.e4/1.e4,ni) 
;    if ni gt 0 then yy(i)=0.
;    cf3(qq).mnvel=cf3(qq).mnvel-yy
;endif
;i=where(cf3.sp2 eq 0 and cf3.jd lt 11413,ni) ;jd 11413 rk33, end of prefix
;if ni gt 0 then begin           ;NO Raw Counts available, use Sine Correction (cf5)
;    cf3(i).mnvel=cf5(i).mnvel
;    cf3(i).errvel = sqrt(cf3(i).errvel^2 + 3.0^2) ;augment errors by 3 m/s in quad 
;    print,'NO Raw Counts found, Sine Correction Applied for: '+cf3(i).obnm
;endif
;endif  ;

;velplot,cf3,'RAW CTS CORRECTED: '+catalog+strupcase(strtrim(strnm,2)),0.6,d1,v1,e1,errcut=2.5,/yrs
velplot,cf3,' '+catalog+strupcase(strtrim(strnm,2)),0.6,d1,v1,e1,errcut=2.5,/yrs

;print,' '
;print,'Saved cf1,cf3, and cf5 in ',vstnm
;print,' '
save,cf1,cf3,file='vstbank/'+'vst'+strlowcase(strtrim(strnm,2))+'.dat'
cfout=cf3
;postscript to laser printer?
if keyword_set(ps) then begin
    if 1-keyword_set(outfile) then outfile = 'jjvank.ps'
    psopen,outfile,xs=8,ys=6,/inch
    if 1-keyword_set(title) then title = catalog+strupcase(strtrim(strnm,2))
    velplot,cf3, title,0.6,/yrs,/nocolor,errcut=2.5  
    psclose
endif
if keyword_set(lazer) then begin
    !p.charsize=2.0
    ps_open,'idl'
    velplot,cf1,catalog+strupcase(strtrim(strnm,2)),0.6,/yrs,/nocolor,errcut=2.5
    ps_close                    ;,/print
    spawn,'lp -d netps1 idl.ps'
;   spawn,'lp -d p433 idl.ps'
    !p.charsize=0
endif



;ans2=''
;read,'Patch the dewar 8 RVs and add prefix velocities? (y/n) ',ans2
;if ans2 eq 'y' then 
if ~keyword_set(nopatch) then patch_d8, star=strlowcase(strnm)
;if ~keyword_set(nopatch) then patch_rawd8, star=strlowcase(strnm), percentile=80
restore,'vstbank/vst'+strlowcase(strnm)+'.dat'
velplot,cf3

return
end


