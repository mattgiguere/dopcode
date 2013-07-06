;PRO DOP_CREATE_IPCF, UPDATE=UPDATE, LIST=LIST
;
; PROCEDURES CALLED: 
; None, but need to edit to put in your paths for the baryfile
;       directory and the vd (files) directory
;
;
; PURPOSE: 
; Reads the vd (files) directory, finds all of the vd's with 
;       good chisq fits and compiles the key information in an 
;       IDL structure.  Used to find the nearest (in time) vd 
;       for an initial wavelength and PSF guess
;
; CALLING SEQUENCE:
; dop_create_ipcf, tag='c', /ctio
;
; INPUTS:
; tag: this is a file tag, not a run tag. It permits the user to 
; carry out several trial Doppler analyses for the same observation
; The file name contains this tag: 'vd'+tag+'iod_rqa06.7453.dat' and 
; 'cd'+tag+'aiod_rqa06.7453' (first pass)  and 
; 'cd'+tag+'biod_rqa06.7453' (second pass)
; 
;
; OPTIONAL INPUTS:
; Update: appends to the existing ipcf.dat file with just a current 
;         observing run (user is queried for the run tag) 
;  
; Written: D. Fischer SFSU Jan 2008
;
;;;;;;;;;;;;;;;;PASS = 1;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dop_create_ipcf, tag=tag, update=update, list=list, goodfit=goodfit, $
                     extra_obs=extra_obs, observatory=observatory,dopenv=dopenv

  if ~keyword_set(observatory) then begin
     observatory=''
     read,'which observatory (ctio, ctio4k, lick, keck): ', observatory
  endif

root_path='/tous/mir7/dop/'
if file_test('/Users/matt/') then root_path='/Users/matt/projects/dopcode/'
if file_test('/Users/debrafischer/') then root_path='/Users/debrafischer/dop1/'
  
     if observatory eq 'ctio4k' then ctio4k=1 else ctio4k=0
     if observatory eq 'ctio' then ctio=1 else ctio=0
     if observatory eq 'keck' then keck=1 else keck=0
     if observatory eq 'lick' then lick=1 else lick=0
     if observatory eq 'het' then het=1 else het=0

  ; Edit the following lines
  if keyword_set(ctio4k) then begin
     vdpath='/tous/mir7/files_df/'
     baryfile = '/tous/mir7/bary/qbcvel.dat'
     ipcf_path=root_path+'ipcf_files/'
     if keyword_set(tag) then ipcf_filename='ipcf_ctio'+tag+'.dat' else $
     	ipcf_filename='ipcf_ctio.dat'
     dewar = 201                ; CTIO designation
     gain = 1.0                 ; CTIO - FITS files have gain applied 
;     gain = 4.5                 ; CTIO - is this correct??? double check 
  endif

  ; Edit the following lines
  if keyword_set(ctio) then begin
     vdpath='/tous/mir7/files/'
     baryfile = '/tous/mir7/bary/qbcvel.dat'
     ipcf_path=root_path+'ipcf_files/'
     ipcf_filename='ipcf_ctio.dat'
     dewar = 200                ; CTIO designation
     gain = 1.0                 ; CTIO - FITS files have gain applied 
;     gain = 4.5                 ; CTIO - is this correct??? double check 
  endif
  if keyword_set(keck) then begin
     vdpath='/tous/mir3/files/'
     baryfile = '/tous/mir3/bary/kbcvel.dat'
     ipcf_path=root_path+'ipcf_files/'
     ipcf_filename='ipcf_keck'+tag+'.dat'
     dewar = 103                ; CTIO designation
     gain = 1.0                ;  Keck fits files 
;     gain = 2.19                ; CTIO - is this correct??? double check 
  endif

  if keyword_set(lick) then begin
     vdpath='/tous/mir1/files/'
     baryfile = '/tous/mir1/bary/bcvel.dat'
     ipcf_path=root_path+'ipcf_files/'
     ipcf_filename='ipcf_lick'+tag+'.dat'
     dewar = 24                ; Lick designation
     gain = 1.264                ;  Lick fits files 
  endif

  if ~keyword_set(goodfit) then goodfit=2.8   ; minimum chisq fit to make it as a reference file

  ; stop editing

  restore,baryfile
  path= vdpath

  if keyword_set(update) then begin
     rtag=' '
     read,'enter the runtag to add: ',rtag
  endif

  print, 'Creating the ipcf.dat file'

;  vdst={obnm:'?', ordt:0, ordob:0, pixt:0, pixob:0, w0:0, wcof:dblarr(dopenv.n_wcof), cts:0L, $
;           z:0.0d, fit:0.0, npix:80, vel:0.0, weight:0.0, $
;           psf:fltarr(121), par:fltarr(npsf)}

  ipcf={obnm:'',iodnm:'',bc:0.0,z:0.0,jd:0.0d,dewar:0,gain:0.0,$
        cts:0L,mnvel:0.0,mdvel:0.0,med_all:0.0,errvel:0.0, mdchi:0.0,$
        nchunk:0,mdpar:fltarr(20),mnpar:fltarr(20),sp1:0.0,sp2:0.0,$
        spst:'',phase:0.0,psfpix:fltarr(17),psfsig:fltarr(17)}

  ipcf_entry = ipcf
  
  if ~keyword_set(list) then begin 
;     iodfiles1=findfile(path+'vd'+tag+'*', count=num1)  &  num2=0
     iodfiles1=findfile(path+'vd'+tag+'iod*', count=num1)
     iodfiles2=findfile(path+'vd'+tag+'HR*', count=num2)
     iodfiles3=findfile(path+'vd'+tag+'MOCKBSTAR*', count=num3)
     if num1 gt 0 and num2 eq 0 then iodfiles=iodfiles1
     if num1 eq 0 and num2 gt 0 then iodfiles=iodfiles2
     if num1 eq 0 and num2 eq 0 then iodfiles=iodfiles3
     if num1 gt 0 and num2 gt 0 and num3 eq 0 then iodfiles=[iodfiles1,iodfiles2]
     if num1 gt 0 and num2 eq 0 and num3 gt 0 then iodfiles=[iodfiles1,iodfiles3]
     if num1 eq 0 and num2 gt 0 and num3 gt 0 then iodfiles=[iodfiles2,iodfiles3]
     if num1 gt 0 and num2 gt 0 and num3 gt 0 then iodfiles=[iodfiles1,iodfiles2,iodfiles3]
     num=n_elements(iodfiles)
  endif ; list not set

  if keyword_set(list) then begin
     readcol,list,iodfiles,f='a20'
     iodfiles=strcompress(path+iodfiles,/remove_all)
     num=n_elements(iodfiles) 
  endif
  if keyword_set(extra_obs) then begin
     iod_extra=strcompress(vdpath+'vd*'+extra_obs,/rem)
     ffex=findfile(iod_extra,count=count_extra)
     if count_extra eq 0 then stop,'not found: '+iod_extra
     if count_extra gt 0 then iodfiles=[iodfiles,ffex[0]]
     num=num+1
  endif

  ipcf_entries = replicate(ipcf, num)

  for i=0,num-1 do begin
     iod_filename = iodfiles[i]
     ;print,iod_filename
     restore, iod_filename    
     ; Check for vdiods that have really bad fits.           
     ; IPCF file.
     nch=n_elements(vd)
     x = where(vd.fit gt goodfit,nx) 
     if (nx gt nch*0.12) then begin
        if nx gt 1 then begin
           print, '*** Not adding ', iod_filename, $
                  ' because it has too many chunks with bad fits!!! ***'
           print, '   You should consider removing this file, since '
           print, '    no one will want to use it.'
        endif
     endif

        nn = strsplit(iod_filename, '/', /extract)
        obn = strsplit(nn[-1],'_',/extract)
        objnm = strmid(obn[0], strlen(tag)+2,strlen(obn[0])-strlen(tag)+2)
        obsnm = obn[1] ; strmid(obn[1],0,strlen(obn[1])-4)
        ipcf_entry.obnm=obsnm
        ipcf_entry.iodnm=objnm

    if strmid(obsnm,0,2) eq 'ac' then tmp_obsnm=strmid(obsnm,1,strlen(obsnm)-1) $
    	else tmp_obsnm=obsnm
	if strmid(obsnm,0,2) eq 'aq' then tmp_obsnm = 'r'+strmid(obsnm,1,strlen(obsnm)-1)
	
     xb=where(bcat.obsnm eq tmp_obsnm,nxb)

     if nxb gt 0 then begin 
        if nxb gt 1 then xb=xb[0]
        ipcf_entry.jd = double(bcat[xb].jd)
        if (ipcf_entry.jd lt 2440000.) then begin
           ipcf_entry.jd = ipcf_entry.jd + 2440000.
        endif
        ipcf_entry.bc=bcat[xb].bc
        ipcf_entry.dewar=dewar
        ipcf_entry.gain=gain
        ipcf_entry.mdchi = median(vd.fit)
        if i eq 0 then ipcf_entries = [ipcf_entry] 
        if i gt 0 then ipcf_entries = [ipcf_entries, ipcf_entry]
     endif

  endfor                        ; each ipcf entry

 ; toss the low chisq fits
  xkeep=where(ipcf_entries.jd gt 0.0 and ipcf_entries.mdchi lt goodfit,nkeep)
  if nkeep eq 0 then stop,'no good fits for ipcf.dat file'
  ipcf_entries=ipcf_entries[xkeep]
                                ; Sort by date.
  ind = sort(ipcf_entries.jd)
  ipcf = ipcf_entries[ind]
    
  num=n_elements(ipcf)

  save, ipcf, file=ipcf_path+ipcf_filename
  if keyword_set(update) then begin
     tmp=ipcf
     restore,ipcf_path+ipcf_filename
     ipcf=[ipcf,tmp]
     ii=sort(ipcf.obnm)
     ipcf=ipcf(ii)
     g=uniq(ipcf.obnm)
     ipcf=ipcf(g)
     jj=sort(ipcf.jd)
     ipcf=ipcf(jj)
     save,ipcf,f=ipcf_path+ipcf_filename
  endif

    print, 'ipcf.dat file saved as: ', ipcf_path+ipcf_filename, $
           ' (', n_elements(ipcf), ' entries)'
    if ~keyword_set(update) then print,$
       'Check to see if the ipcf.dat file needs to be copied to a new name'
;stop    
end  ;pro

