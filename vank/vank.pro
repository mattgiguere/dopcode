pro vank,strnm,tmptp,lbl,cf1,cf3,mct=mct,lazer=lazer,maxchi=maxchi,$
noclean=noclean,cf=cfout,ps=ps,title=title,dfd=dfd,order=order,$
outfile=outfile,vdarr=vdarr, reject=reject,dircf=dircf,fit_key=fit_key,$
doptest=doptest,nopatch=nopatch, keck=keck,ctio=ctio,lick=lick, test=test

; vank,'moon','c','vdcnso_',cf1,cf3,mct=8000.  
; lbl is everything before the starnm in the vd files
; vank, '128620','adf','vdadf', cf1,cf3,mct=10000,/ctio, fit_key=1

;fit_key = 1 for stars with CHIRON templates
;fit_key = 1 for stars with lick templates
;fit_key = 2 for stars with keck templates

if ~keyword_set(tmptp) then tmptp='b'
if ~keyword_set(lbl) then lbl='vdbnso'
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

;mk_cf,strnm,tag=tmptp
  m7_psfpix = [0.0, -4.0, -3.6, -3.0, -2.4, -1.8, -1.2, -0.8, -0.4, 0.4, 0.8,  1.2,  1.8,  2.4, 3.0, 3.6, 4.0]
  m7_psfsig = [1.05, 0.0,  0.8,  0.7,  0.7,  0.7,  0.7,  0.6,  0.5, 0.5, 0.6,  0.7,  0.7,  0.7, 0.7, 0.8, 0.0]

;   m7_psfpix = 1.2*[0.00, -2.2,-1.7, -1.2, -0.7,-0.3, 0.3, 0.7, 1.2, 1.7, 2.2]
;   m7_psfsig = 1.2*[1.10,  0.6, 0.6,  0.6,  0.6, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6]

vstdsk = '/Users/debrafischer/dop2/vank/vstbank/'
if keyword_set(ctio) then begin
   restore,'/tous/mir7/bary/qbcvel.dat'
   if ~keyword_set(dfd) then dfd='/tous/mir7/files_df/'
   if keyword_set(test_dir) then dfd=dfd+test+'/'
   fildsk = '/tous/mir7/files_df/'+lbl+strnm+'_'
   cfdsk = '/tous/mir7/cf/'
endif

if keyword_set(lick) then begin
   restore,'/mir1/bary/bcvel.dat'
   if ~keyword_set(dfd) then dfd='/mir1/files/'
   fildsk = '/mir1/files/'+lbl+strnm+'_'
   cfdsk = '/mir1/cf/'
endif

if keyword_set(keck) then begin
   restore,'/mir3/bary/kbcvel.dat'
   if ~keyword_set(dfd) then dfd='/mir3/files/'
   fildsk = '/mir3/files/'+lbl+strnm+'_'
   cfdsk = '/mir3/cf/'
endif
grnm=dfd+lbl+strnm
ff=file_search(grnm+'*',count=count)
;x=where(strmid(ff,strlen(ff)-5,5) ne '_orig')
;ff=ff[x]

tmp={cfstr,obnm:'?',iodnm:'?',bc:0.,z:0.,jd:double(0.),dewar:50,gain:1.3, $
	 cts:long(0),mnvel:0.,mdvel:0.,med_all:0.,errvel:0.,mdchi:0.,nchunk:0, $
	 mdpar:fltarr(20),mnpar:fltarr(20),sp1:0.,sp2:0.,spst:'?',phase:0., $
	 psfpix:[m7_psfpix[0:10], 0.00, 0.00, 0.00,0.00],$
	 psfsig:[m7_psfsig[0:10], 0.00, 0.00, 0.00,0.00]}
	cf=replicate(tmp,count)

for i = 0, count-1 do begin
   x1=strpos(ff[i],'_',/reverse_search)
   obnm=strmid(ff[i],x1+1,strlen(ff[i])-x1)
   if strmid(obnm,0,1) eq 'a' then x=where(bcat.obsnm eq strmid(obnm,1,strlen(obnm)-1),nx) $
   else x=where(bcat.obsnm eq obnm,nx)
   if nx eq 0 then stop,' no matches in bcat'
   if nx gt 1 then stop,' more than one match in bcat'
   cf[i].obnm=obnm
   cf[i].bc=bcat[x].bc
   cf[i].jd=bcat[x].jd
endfor

if strnm eq '10700' then begin 
	ss=where(cf.jd lt 16109.0 or cf.jd gt 16110.0,nss) ;exclude 120630
	cf=cf[ss]
;	print,cf.obnm
;	stop
endif

save,cf,f=cfdsk+strupcase(strtrim(strnm,2))+'_'+tmptp+'.dat'
cfnm    = cfdsk  + strupcase(strtrim(strnm,2))+'_'+tmptp+'.dat'

if keyword_set(dircf) then cfnm = dircf
vstnm   = vstdsk + 'vst'+strlowcase(strtrim(strnm,2))+'.dat'
vnm     = strnm

if n_elements(mct) eq 1 then if mct gt 0 then begin
    print,'Using input value for minimum acceptable photons in exposure:  '+strtrim(mct,2)
    mncts=mct
endif

ff_cf=findfile(cfnm,count=cfcount)

if cfcount gt 0 then restore,cfnm else message,'CF file '+cfnm+' not found'

cff=cf

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
        ordr=[1, 12]
        pixr=[50, 3970]
    end
    tape eq 'rq' : begin
;        ordr=[5, 20]
;        pixr=[280,1640]
        ordr=[12,29]
        pixr=[80,3040]
    end
    tape eq 'ac' : begin
		if keyword_set(order) then ordr=order else ordr=[13,32]
		pixr=[80, 3040]
    end
    tape eq 'ch' : begin
		if keyword_set(order) then ordr=order else ordr=[13,32]
		pixr=[80, 3040]
    end
    tape eq 'em' and cff[0].jd ge 2.44d6+13480: begin
        ordr=[0, 14]
        pixr=[50, 3890]
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
        ;dumfit=3.*median([vd0.fit,vd1.fit])
        ;dumwt=0.5*median(vd0.weight)
        ;ind=where(vd0.fit lt dumfit and vd1.fit lt dumfit and $
        ;          vd0.weight gt dumwt and $
        ;          vd0.pixt gt pixr(0) and vd0.pixt lt pixr(1) and $
        ;          vd0.ordt ge ordr(0) and vd0.ordt le ordr(1),nchunk)
		ind=where(vd0.pixt gt pixr[0] and vd0.pixt lt pixr[1] and $
				  vd0.ordt ge ordr[0] and vd0.ordt le ordr[1], nchunk)
        cf1(0).mnvel=mean(vd0(ind).vel)  &  cf1(0).mdvel=median(vd0(ind).vel)
        cf1(1).mnvel=mean(vd1(ind).vel)  &  cf1(1).mdvel=median(vd1(ind).vel)
        cf1(0).errvel=stdev(vd0(ind).vel)/sqrt(nchunk-1)
        cf1(1).errvel=stdev(vd1(ind).vel)/sqrt(nchunk-1)
        cf1(0).mdchi=median(vd0(ind).fit)
        cf1(1).mdchi=median(vd1(ind).fit)
        cf1(0).nchunk=nchunk
        cf1(1).nchunk=nchunk
        cf1(0).cts=median(vd0(ind).cts)
        cf1(1).cts=median(vd1(ind).cts)
        cf1.jd=cf1.jd-2440000.
   endif else begin
   		if keyword_set(fit_key) then $
   			vel,vnm,cff,lbl,cf1,vdarr,ordr=ordr,pixr=pixr,$
   				upvel=1,dwr=dwr, vdpath=dfd,mincts=mncts,maxchi=maxchi,ifit=fit_key, $
   				noclean=noclean,doptest=doptest else $
			vel,vnm,cff,lbl,cf1,vdarr,ordr=ordr,pixr=pixr,mincts=mncts,maxchi=maxchi,dwr=dwr,$
	    		vdpath=dfd,ifit=ifit,noclean=noclean,doptest=doptest
	end ;else
	endif else begin
    	print,'Less than two observations'
    return
endelse

bv=0.                  
cf3=cf1
velplot,cf3,0.006,d1,v1,e1,errcut=2.5

if keyword_set(order) then sufx='_ordr' else sufx=''

vstnm1='vst'+strlowcase(strtrim(strnm,2))+sufx+'.dat'
print,'vstnm1: ',vstnm1
save,cf1,cf3,file='vstbank/'+vstnm1
ck_shutter,star=strnm, vstnm=vstnm1   ;reject observations with midpoint time errors

cfout=cf3

;postscript to laser printer?
if keyword_set(ps) then begin
    if 1-keyword_set(outfile) then outfile = 'vank.ps'
    psopen,outfile,xs=8,ys=6,/inch
    if 1-keyword_set(title) then title = strupcase(strtrim(strnm,2))
    velplot,cf3, 0.006,/yrs,/nocolor,errcut=2.5  
    psclose
endif
if keyword_set(lazer) then begin
    !p.charsize=2.0
    ps_open,'idl'
    velplot,cf1,0.006,/yrs,/nocolor,errcut=2.5
    ps_close                    ;,/print
    spawn,'lp -d netps1 idl.ps'
;   spawn,'lp -d p433 idl.ps'
    !p.charsize=0
endif

restore,'vstbank/'+vstnm1
velplot,cf3;,/fitline

return
end


