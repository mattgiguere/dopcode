pro vel, starname, cf, label, outcf, vdarr, med_all, acof,             $
           plot=plot, pixr=pixr, ordr=ordr, mincts=mincts,dwr=dwr,     $
	   vdpath=vdpath, ifit=ifit, maxchi = maxchi,trend=trend,      $ 
	   detrend=detrend, noclean=noclean,newvd=newvd,weight=weight, $ 
	   post_fix=post_fix,upvel=upvel,sp2=sp2,noplot=noplot,doptest=doptest
;
;  Reads all "VD" structures having the label, "LABEL" (i.e., 'vdE') on disk 
;  for the star, "STARNAME", as specified in file, "CF",
;  which must be restored prior to operation.
;  It creates an output velocity structures, VDCUBE and VELST, of 
;  times of observation and average velocity for each observation. 
;
;INPUT:
;  starname (input string)    examples: '509' or '4983' or 'GL897'
;  cf       (input structure) examples: restore,'/mir1/ew/doppler/cf.GL380'
;  label    (string)          VD label  (= 'vd', 'vdA', 'vdB', ...)
;
;OPTIONAL OUTPUT:
;  vdarr     (output structure)  contains JD, avg. velocities, etc.
;  cfout     (output structure)  trimmed cf (excluding bad obs., etc.)
;  med_all   (float)             The median velocity of all chunks
;
;OPTIONAL KEYWORDS:
;  plot        (keyword)       /plot produces Velocity vs. Time .
;  pixr        ( intarr(2) )   Pixel Range -   default: pixr=[160,640] 
;  ordr        ( intarr(2) )   Order range -   default: ordr=[6,21]
;  mincts      ( long )        Minimun acceptable counts
;  dwr         (int)           Dewars to exclude
;  vdpath      (string)        directory path for VD files (/mir1/files/)
;  ifit        ( keyword )     (= 1 or 2: vd.ifit,  vd.sfit for Chi-Sq values)
;                              Use ifit=1 when last pass employed /PSF.
;                              Use ifit=2 when last pass employed /ZPASS
;  maxchi      ( keyword )     toss observations having median(chi) > maxchi
;  newvd        (keyword)      set to 1 for newvd (for asymmetry detrend)
;  upvel        (keyword)      set to 1 to update velocities from "vd.z" and "cf.bc"
;ROUTINES USED:
;   VDCUBE    to read in VD structure of each observation.
;   VDCLEAN   to purge VDCUBE of chunk sets having poor Chisq and Slopeweight.
;   DETREND   remove correlation of velocity vs. IP params
;   VELPLOT   to plot Velocity vs. Time
;
;HISTORY:
;   Spinoff from VEL.PRO  circa May 8 - 11, 1993  R.P.B.
;   Modified  Apr 17, 1994 RPB,GWM for CF driven version.
;   
;Modified for VDCUBE.PRO  May 1994 for VDARR output. GWM
;
IF n_params () le 2 then begin
  print,'-------------------------------------------------------------------'
  print,' SYNTAX:'
  print,' '
  print,' IDL> restore,''/mir1/ew/doppler/cf.GL380''    ;get CF structure'
  print,' '
  print,' IDL> vel,star, cf, label, [cfout, vdarr, cfout, med_all, acof],    '
  print,'            plot=plot, pixr=[a,b], ordr=[a,b], mincts=mincts, dwr=dwr '
  print,'           vdpath=''vdpath/'', ifit=ifit, maxchi=maxchi, '
  print,'           noclean=noclean, newvd=newvd'
  print,' '
  print,'Example: vel,''8086'',cf,''vdE'',/plot,dwr=13'
  print,'-----------------------------------------------------------------'
  RETURN
ENDIF
;
;SET CONSTANTS and DEFAULTS
;    starname = strupcase(strtrim(starname,2))       ;trim blanks from starname)
    starname = strtrim(starname,2)                  ;trim blanks from starname)
    c = 2.997925d8                                  ;speed of them photons
;    minwt = 0.6
  if n_elements(upvel) eq 1 then if upvel eq 1 then $
     print,'Updating velocites with new barycentric corrections.'
  if not keyword_set(pixr)    then pixr = [000,4000]   ;pixel range
  if not keyword_set(ordr)    then ordr = [0,60]      ;order range
  if not keyword_set(mincts)  then mincts = 1         ;min. required counts
  if not keyword_set(dwr)     then dwr = 0            ;reject dewars
  if not keyword_set(ifit)    then ifit = 0           ;0 (no /psf) 1 (/psf)
  if not keyword_set(maxchi)  then maxchi = 3         ;default chi-sq value
  if not keyword_set(newvd)   then newvd = 0          ;default to old vd
  if keyword_set(vdpath) then begin
     vdpath=strtrim(vdpath,2)                         ;Trim VD dir path
     len = strlen(vdpath)                             ;check last char...
     lastchar = strmid(vdpath,len-1,1)                ;It should be '/'
     if lastchar ne '/' then vdpath = vdpath + '/'    ;put "/" at end
  end
  label=[label]
  dwr = [dwr]                                    ;rejected dewars
  mincts = long(mincts)                          ;required cts

inpcf=cf

;VDCUBE: get VDARR(obs#,chunk#) STACK OF VD PANCAKES.
  vdcube,starname, inpcf, label, vdarr, outcf, nchunk, med_all,  $
         pixr=pixr, ordr=ordr, mincts=mincts, dwr=dwr, ifit=ifit,    $
         vdpath=vdpath, maxchi=maxchi,weight=weight,upvel=upvel,sp2=sp2
;
   tagnam = tag_names(vdarr)
   fitdex = first_el(where(tagnam eq 'FIT')) ;Def. Tag_name index of chi-!
   veldex = first_el(where(tagnam eq 'VEL')) ;Def. Tag_name index of chi-!
   if ifit eq 1 then begin
	fitdex=first_el(where(tagnam eq 'FIT'))
;        veldex = first_el(where(tagnam eq 'IVEL')) ;Def. Tag_name index of chi-!
        veldex = first_el(where(tagnam eq 'VEL')) ;Def. Tag_name index of chi-!
	if keyword_set(doptest) then veldex = first_el(where(tagnam eq 'VEL'))
   endif

  if ifit eq 2 then begin
	fitdex=first_el(where(tagnam eq 'SFIT'))
        veldex = first_el(where(tagnam eq 'SVEL')) ;Def. Tag_name index of chi-!
   endif

  if ifit eq 3 then begin
	fitdex=first_el(where(tagnam eq 'FIT1'))
        veldex = first_el(where(tagnam eq 'VEL1')) ;Def. Tag_name index of chi-!
   endif

   if n_elements(sp2) eq 1 then if sp2 eq 1 then veldex = first_el(where(tagnam eq 'SP2'))
;
   orddex = first_el(where(tagnam eq 'ORDER')) ;Def. Tag_name index of ORDER!
   if orddex eq -1 then orddex = first_el(where(tagnam eq 'ORDT'))
;

  if not keyword_set(noclean) then begin
    if not keyword_set(noplot) then VDCLEAN,vdarr,ifit,nchunk,sp2=sp2    ;Reject chunks w/ poor chisq and slopeweight
    if keyword_set(noplot) then VDCLEAN,vdarr,ifit,nchunk,sp2=sp2,/noplot    ;Reject chunks w/ poor chisq and slopeweight
;    print,'Chunk Set zero-pts NOT SET TO MEDIAN in VDCLEAN (/nozero).'
  end

  numobs = n_elements(vdarr(*,0))
   dum  = outcf.obnm
   if min(outcf.jd) gt 2440000 then dum=outcf.JD-2440000. else dum=outcf.JD
   outcf.JD    = dum
   outcf.MED_ALL= med_all	
;

;  OBSERVATIONS LOOP
   FOR ob = 0,numobs-1 do begin                  ;LOOP THRU OBSERVATIONS
;    OUTLIER REJECTION: reject chunks with far out Velocities
     velo = reform(vdarr(ob,*).(veldex))
     resid = abs(velo - median(velo)) ;residual velocity from median
     indres = sort(resid)             ;Sort residuals
     numind = n_elements(resid)
     thr = 0.99
     gd1 = where (vdarr(ob,*).weight lt 100 and vdarr(ob,*).weight gt 0.00)
     gd2 = indres(0:thr*(numind-1))    ;indices of chunks: exclude worst 1% of vel residuals
     gd=[gd1,gd2]
     i=sort(gd)
     gd=gd(i)
     g=uniq(gd)
     gd=gd(g)	
;
     velo = vdarr(ob,gd).(veldex)                     ;Vels in this obs, having low chi
     wt =   vdarr(ob,gd).weight
     fit =  vdarr(ob,gd).(fitdex)
     vind = where(fit lt 99) 
     velo=velo(vind)
     wt  =  wt(vind)
     outcf(ob).mdvel = median(velo)               ;median vel
     outcf(ob).mnvel = total(velo*wt)/total(wt)   ;weighted mean vel
     outcf(ob).MDCHI = median(vdarr(ob,gd).(fitdex))     ;median chi-sq
     outcf(ob).CTS   = median(vdarr(ob,gd).cts)          ;median cts
     outcf(ob).NCHUNK = n_elements(velo)
     outcf(ob).z = median(vdarr(ob,gd).z)        ;median z jun 2006
     outcf(ob).mnvel = outcf(ob).mnvel+ outcf(ob).z*outcf(ob).bc ;new BC jun 2006
     
;
;     for i=0,14 do begin
; comment the next 4 lines out for vank bombs... fischer 11-98

if starname eq '4112' or starname eq '5011' or $
     starname eq '5553' or starname eq '5933' then goto, jumpover
;     for i=0,n_elements(vdarr(0,0).iparam)-1 do begin
;       outcf(ob).mdpar(i) = median(vdarr(ob,gd).iparam(i))
;      outcf(ob).mnpar(i) = mean(vdarr(ob,gd).iparam(i))
;     end
jumpover: 
;
;    ERROR in the MEAN (see Butler et al. 1996 PASP, eqn.4)
;    WT  must be 1/sigma_vel^2   having units of 1/(m/s)^2
     outcf(ob).errvel = 1./sqrt( total(wt) )     ;weighted error in mean    
   ENDFOR                                        ;end observation loop
;

;   if keyword_set(detrend) then begin
;     DETREND,vdarr,outcf,acof,plot=plot,newvd=newvd
;   end
   print,' '
   print,'Number of Final Observations Analyzed:',numobs
   print,' '
;

;PLOTTING SECTION
  if keyword_set(plot) then begin 
    ;VELPLOT,vst,starname,/yrs,
;    vst.mnvel = vst.mnvel - median(vst.mnvel)
    velplot, outcf,starname,0.6,yra=[-50,50],/yrs,/lzr ;,/nobox ;,/nocolor
end
return
end
