pro vdcube, starname, inpcf, label, vdarr, outcf, nchunk, med_all,   $
           pixr=pixr, ordr=ordr, mincts=mincts, dwr=dwr,        $
           vdpath=vdpath, ifit=ifit, maxchi = maxchi,           $
           nozero=nozero,weight=weight,upvel=upvel,sp2=sp2,     $
           ord_skip=ord_skip,gpix=gpix                                      
;
;  Reads all "VD" structures having the label, "LABEL" (i.e., 'vdE') on disk 
;  for the star, "STARNAME", as specified in structure, "INPCF",
;  which must be RESTOREd prior to operation, or generated with CF.PRO.
;  It creates an output VD cube, VDARR, containing the VD's which are found and
;  which meet the criteria of mincts, dwr, and maxchi.  VDARR(3,48).iparam(0)
;  is the central Gaussian width for observation 3, chunk 48.
;
;INPUT:
;  starname (input string)    examples: '509' or '4983' or 'GL897'
;  inpcf    (input structure) examples: restore,'/mir1/ew/doppler/cf.GL380'
;  label    (string)          VD label  (= 'vd', 'vdA', 'vdB', ...)
;
;OUTPUT:
;  vdarr    (output structure)  contains an array of VD structures, for 
;                               each observation and chunk within it.
;  outcf    (output structure)  Extracted CF of found and useful observations
;  nchunk   (integer)           # of chunks in each vd structure per obs.
;
;OPTIONAL:
;  pixr        ( intarr(2) )   Pixel Range -   default: pixr=[180,620] 
;  ordr        ( intarr(2) )   Order range -   default: ordr=[6,21]
;  mincts      ( long )        Minimun acceptable counts
;  dwr         (int)           Dewars to exclude
;  vdpath      (string)        directory path for VD files (/mir1/files/)
;  ifit        ( keyword )     (= 1 or 2: vd.ifit,  vd.sfit for Chi-Sq values)
;  maxchi      ( keyword )     toss observations having median(chi) > maxchi
;                              Use ifit=1 when last pass employed /PSF.
;                              Use ifit=2 when last pass employed /ZPASS
; nozero (keyword)             Don't add constant to each chunk-set ==> zero avg. vel.
;                              Default: subtract median of chunk set.
; weight (keyword)             Replace weights in "VD"
; upvel  (keyword)             Recalculate velocities using "VD.Z" and "CF.BC"
;
;Spinoff from VEL.PRO  circa May 8 - 11, 1993  R.P.B.
;Modified  Apr 17, 1994 RPB,GWM for CF driven version.
;Modified as VDCUBE.PRO  May 1994 for VDARR output.
;
IF n_params () lt 4 then begin
  print,'-------------------------------------------------------------------'
  print,' SYNTAX:'
  print,' '
  print,' IDL> restore,"/mir1/ew/doppler/cf.GL380"'    ;get CF structure
  print,' '
  print,' IDL> vdcube,star, inpcf, label, vdarr, outcf'
  print,'            pixr=[a,b], ordr=[a,b], mincts=mincts, dwr=dwr '
  print,'            vdpath=''vdpath/'', ifit=ifit, maxchi=maxchi, '
  print,' '
  print,'      where ''star'' is e.g., ''509'' or ''GL380'' '
  print,'-----------------------------------------------------------------'
  RETURN
ENDIF
;   
;Yale
   vdpth0 = '/tous/mir7/files_df/'                    
   vdpth1 = '/tous/mir7/files_df/'
   vdpth2 = '/tous/mir3/files/'

;SET CONSTANTS and DEFAULTS
;    starname = strupcase(strtrim(starname,2))       ;trim blanks from starname)
    starname = strtrim(starname,2)                   ;trim blanks from starname)
    cf=inpcf                                         ;internal var, CF
    c = 2.99792458d8                                 ;speed o' dem der photons
  if not keyword_set(pixr)     then pixr = [0,2000]  ;pixel range
  if not keyword_set(ordr)     then ordr = [40,53]   ;order range
  if not keyword_set(ord_skip) then ord_skip = [-1]  ;toss orders
  if not keyword_set(mincts)   then mincts = 0.1     ;min. required counts
  if not keyword_set(dwr)      then dwr = 0          ;reject dewars
  if not keyword_set(ifit)     then ifit = 0         ;default chi-sq value
  if not keyword_set(maxchi)   then maxchi = 100     ;default chi-sq value
  if keyword_set(vdpath) then begin
     vdpath=strtrim(vdpath,2)                       ;Trim VD dir path
     len = strlen(vdpath)                           ;check last char...
     lastchar = strmid(vdpath,len-1,1)              ;It should be '/'
     if lastchar ne '/' then vdpath = vdpath + '/'  ;put "/" at end
  end


  label = strtrim(label,2)                        ;Trim VD label
  vdstar  = label + starname
  dwr = [dwr]                                    ;rejected dewars
  mincts = long(mincts)                          ;required cts
  highcts = 100000                               ;high counts on Hamilton
  if strmid(cf(0).obnm,0,2) eq 'rk' then highcts=1000000 ;high cts on Keck
  if strmid(cf(0).obnm,0,2) eq 'rj' then highcts=1000000 ;high cts on Keck
  if strmid(cf(0).obnm,0,2) eq 'em' then highcts=4d3^2 ;high cts on Keck
  if strmid(cf(0).obnm,0,2) eq 'ru' then highcts=300000  ;high cts on AAT
  if strmid(cf(0).obnm,0,2) eq 'rv' then highcts=1000000 ;high cts on UVES
;
;

;DEWAR TOSS                                      ;to see who kicks off
	IF dwr(0) ne 0 then begin                ;begin dewar toss
          for qq = 0,(n_elements(dwr)-1) do begin
              ind = where(cf.dewar ne dwr(qq))
              cf  = cf(ind)
              print,'Tossed observations made with dewar ',dwr(qq)
          endfor
          print,'in loop '
	ENDIF
        numob = n_elements(cf)
        qq = -1   ;counter to find first "VD"

;;;path=vdpth0
;;for qq=0,numob-1 do vdname = path + vdstar + '_'+ cf(qq).obnm          ;VD disk name

   repeat begin
        qq=qq+1
;print,qq

;RESTORE FIRST VD FOR TAG_NAMES (= 0th VD from CF)
        IF keyword_set(vdpath) then path = vdpath    $ 
             ELSE begin                                 ;if VDPATH not specified 
               rx = strmid(cf(qq).obnm,0,2)              ; one of:  rh ra rk ...
;print,rx
               CASE rx OF
                  'ra': path = vdpth0
                  'rb': path = vdpth0
                  'rc': path = vdpth1
                  'rd': path = vdpth0
                  're': path = vdpth0
                  'rf': path = vdpth0
                  'rg': path = vdpth0
                  'rh': path = vdpth0
                  'ri': path = vdpth0
                  'rj': path = vdpth2
                  'rk': path = vdpth2
                  'rm': path = '/mir5/files/'
                  'rq': path = '/tous/mir7/files/'
                  'ch': path = '/tous/mir7/files_df/'
                  'ac': path = '/tous/mir7/files_df/'
                  'rs': path = vdpth0
                  'rt': path = vdpth0
                  'rx': path = vdpth0
                  'rz': path = vdpth0
                  'em': path = '/mir3/files/'
               ELSE: print,rx+' VDCUBE: Directory Path to VD was not found!'
               ENDCASE
             ENDELSE

	vdname = path + vdstar + '_'+ cf(qq).obnm          ;VD disk name

   endrep until (first_el(findfile(vdname)) eq vdname)     ;Must find first VD
          if first_el(findfile(vdname)) ne vdname then begin   ;File doesn't exist
            print,'VDCUBE: Did not find first VD file: ',vdname
            return
          endif

        restore,vdname
        vd[129].weight=0.0

	tagnam = tag_names(vd)
        veldex = first_el(where(tagnam eq 'VEL'))
        if ifit eq 1 then begin
             veldex = first_el(where(tagnam eq 'VEL'))
             if veldex lt 0 then veldex = first_el(where(tagnam eq 'VEL'))
        endif
        if ifit eq 2 then begin
             veldex = first_el(where(tagnam eq 'SVEL'))
             if veldex lt 0 then veldex = first_el(where(tagnam eq 'VEL'))
        endif
        if ifit eq 3 then begin
             veldex = first_el(where(tagnam eq 'VEL1'))
             if veldex lt 0 then veldex = first_el(where(tagnam eq 'VEL'))
        endif
;	if n_elements(sp2) eq 1 then if sp2 eq 1 then veldex = first_el(where(tagnam eq 'SP2'))

	if n_elements(upvel) eq 1 then if upvel eq 1 then $
	   vd.(veldex)=double(vd.z)*c+cf(qq).bc

	tagnam=tag_names(vd)
        ordindex=first_el(where(tagnam eq 'ORDER'))
        if ordindex lt 0 then ordindex=first_el(where(tagnam eq 'ORDOB'))
;PB Kludge to patch in modern GM weights
	if n_elements(weight) eq n_elements(vd) then vd.weight=weight
;
;VD INDICES USED: indices satisfying Pixel and Order criteria
      WPOind = where(vd.pixt  ge pixr[0]-.1 and $   
                     vd.pixt  le pixr[1]+.1 and $
                     vd.(ordindex) ge ordr[0]-.1 and $
                     vd.(ordindex) le ordr[1]+.1, nchunk )
      vd = vd[WPOind]
;PB Kludge to toss individual orders, typically tellurics near 6000 angstroms  23Mar2003
      if ord_skip(0) gt -1 then $
          for zz=0,n_elements(ord_skip)-1 do vd=vd(where(vd.(ordindex) ne ord_skip(zz)))
      nchunk=n_elements(vd)
      print,' '
      ord1 = strtrim(string(min(vd.(ordindex))),2)
      ord2 = strtrim(string(max(vd.(ordindex))),2)
      pix1 = strtrim(string(min(vd.pixt)),2)
      pix2 = strtrim(string(max(vd.pixt)),2)
      print,'                      Order Range: ',ord1,'-',ord2
      print,'                      Pixel Range: ',pix1,'-',pix2
      print,'                      # of Chunks: ',strtrim(string(fix(nchunk)),2)
      initvel = median(vd.(veldex))              ;ruf initial velocity      
;
;
;VDARR = replicate({vd},numob,nchunk)    ;All VD's:  VD=vdarr(ob,chunk)
;TAG_NAMES ESTABLISHED
      tagnam = tag_names(vd)
      fitdex = first_el(where(tagnam eq 'FIT'))    ;Def. Tag_name index of chi-sq
      if ifit eq 1 then fitdex=first_el(where(tagnam eq 'FIT'))
      if ifit eq 2 then fitdex=first_el(where(tagnam eq 'SFIT'))
      if ifit eq 3 then fitdex=first_el(where(tagnam eq 'FIT1'))
VDARR = replicate(vd[0],numob,nchunk)    ;All VD's:  VD=vdarr(ob,chunk)
;
;PRINT HEADER
     print,' '
     print,'--------------------------------------------------------------------------
     print,' Restored File  Median Velocity     Median Photons  Median'
     print,'                 *All Chunks*        *All Chunks*    Fit'
     print,'-------------------------------------------------------------------------'
;
;LOOP THROUGH ALL OBSERVATIONS
     numgood = 0                    ;# of good observations
     gdind = [-1]                   ;index of good observations
     velarr = fltarr(numob)         ;array of median velocity
     FOR ob = qq, numob-1 do begin
        IF keyword_set(vdpath) then path = vdpath $
             ELSE begin                             ;VDPATH not specified 
               rx = strmid(cf(ob).obnm,0,2)          ; one of:  rh ra rk ...
               CASE rx OF
                  'ra': path = vdpth0
                  'rb': path = vdpth0
                  'rc': path = vdpth1
                  'rd': path = vdpth0
                  're': path = vdpth0
                  'rf': path = vdpth0
                  'rg': path = vdpth0
                  'rh': path = vdpth0
                  'ri': path = vdpth0
                  'rj': path = vdpth2
                  'rk': path = vdpth2
                  'rm': path = '/mir5/files/'
                  'rq': path = '/mir7/files/'
                  'ch': path = '/tous/mir7/files/'
                  'ac': path = '/tous/mir7/files/'
                  'rs': path = vdpth0
                  'rt': path = vdpth0
                  'rx': path = vdpth0
                  'rz': path = vdpth0
                  'em': path = vdpth2
;                 path = paths(1, where( paths(0,*) eq rx ) )
               ENDCASE
        END
        vdname = path + vdstar + '_'+ cf[ob].obnm

        dum = first_el(findfile(vdname))
        IF dum ne vdname then begin
           messag1 = 'NOT FOUND.'
           messag2 = ' '
           messag3 = ' '
           goto, JUMP1                      ;go to bottom of loop
        ENDIF
;
        restore,vdname
        vd[129].weight=0.0

	if n_elements(upvel) eq 1 then if upvel eq 1 then $
	   vd.(veldex)=double(vd.z)*c+cf(ob).bc

;PB Kludge to patch in modern GM weights
	if n_elements(weight) eq n_elements(vd) then vd.weight=weight
        vd = vd(WPOind)                       ;VD CHUNKS (Desired Wt, Pix, Ord) 
;PB Kludge to toss individual orders, typically tellurics near 6000 angstroms  23Mar2003
        if ord_skip(0) gt -1 then $
            for zz=0,n_elements(ord_skip)-1 do vd=vd(where(vd.(ordindex) ne ord_skip(zz)))
        medvel = median(vd.(veldex))          ;quick-look results
        medcts = long(median(vd.cts)) 
        medchi = median(vd.(fitdex))
        velarr(ob) = medvel

;
;       LOW-COUNT OBSERVATIONS TOSSED    (only if  /mincts keyword set)
        IF medcts lt mincts then begin
           messag1 = 'TOSSED Due'
           messag2 = 'to Low Counts:  '+strtrim(medcts,2)
           messag3 = ' '
           goto, JUMP1                       ;go to bottom of loop
        ENDIF
;
;       HIGH-COUNT OBSERVATIONS TOSSED
        IF medcts gt highcts then begin
           messag1 = 'TOSSED Due'
           messag2 = 'to High Counts:  '+strtrim(medcts,2)
           messag3 = '  '
           goto, JUMP1                       ;go to bottom of loop
        ENDIF

;       HIGH CHI-SQ OBSERVATIONS TOSSED  (only if /maxchi keyword set)
        IF medchi gt maxchi then begin
           messag1 = 'Tossed Due'
           messag2 = 'High Chi-Sq: '+strtrim(medchi,2)
           messag3 = '  '
           goto, JUMP1                       ;go to bottom of loop
        ENDIF

        IF medchi eq 0 then begin
           messag1 = 'Tossed Due'
           messag2 = '0 Chi-Sq'
           messag3 = '  '
           goto, JUMP1                       ;go to bottom of loop
        ENDIF
;
;       Execute here if VD was both found AND met chisq, mincts criteria.
        numgood = numgood+1                  ;cumul # of good obs (1 base)
        gdind = [gdind,ob]                    ;index of good obs
        ;messag1 = strtrim(fix(medvel-initvel),2)
        messag1 = string(medvel-initvel) ; dave
        messag2 = strtrim(medcts,2)
        messag3 = strmid(strtrim(medchi,2),0,5)
        ntg=n_tags(vd)
;        if ntg gt 14 then vd = del_tag2(vd,'par')
        vdarr(numgood-1,*) = vd         ;stuff VD into VDARR
 JUMP1: format = '(A18,3x,A11,A22,A10)'
        pos = strpos(vdname,'vd')
;        vdtrim = strmid(vdname,pos,17)
        vdtrim = strtrim(cf(ob).obnm,2)
        print,format=format,vdtrim,messag1,messag2,messag3
     ENDFOR ;       End Loop thru all observations ( ob = 0,numob-1)

;
     print,'---------------------------------------------------------'
;
     if numgood le 1 then begin
        print,'N O   O B S E R V A T I O N S   W E R E   R E T A I N E D!'
        return
     endif
     gdind = gdind(1:numgood)                ;trim -1 from 1st element
     outcf = cf(gdind)                         ;extract "good" obs from CF
     vdarr = vdarr(0:numgood-1,*)              ;extract "good" vdarr's
     velarr = velarr(gdind)
     avvel = mean(velarr)
     sig = stddev(velarr)
     format = '(A30,3x,F12.1,A5)'
     print,format=format,'    AVERAGE          ',avvel,' m/s' 
     print,format=format,'      SIGMA          ',sig,' m/s'
     print,' '
;
     med_all = median(vdarr(*,*).(veldex))      ;median velocity of all chunks
     
return
end
