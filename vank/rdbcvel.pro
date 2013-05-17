pro rdbcvel,starname,cf,obdsk,logfile=logfile,noprint=noprint, obdsk=obdsk1,obnm_filter=obnm_filter,dewar=dewar,svlist=svlist
if keyword_set(obdsk1) then obdsk = obdsk1
;This code reads bcvel and constructs a cf data structure with the
; necessary information
;
;starname (input string)    examples: '509' or '4983' or 'GL897'
;logfile (keyword string)   This allows on-line observation log sheets
;			      other that the default (bcvel.ascii)
;			      to be searched.
;obnm_filter (input list of strings) Optional list of particular observation names we're interested in. 
;                                    All observations for this star not in this list are ignored.
;dewar                       Specifies to only load observations taken with this dewar.

;Created Aug 5, 1993  R.P.B.
;
               print, '*** *** dave. this is a kludge n rdbcvel. take it out once done *** ***' ; 3/17/06
               obd0 = '/data/ctio/iodspec/'

   m7_psfpix = 1.2*[0.00, -2.2,-1.7, -1.2, -0.7,-0.3, 0.3, 0.7, 1.2, 1.7, 2.2]
   m7_psfsig = 1.2*[1.10,  0.6, 0.6,  0.6,  0.6, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6]

if ~keyword_set(logfile) then begin
       starname=strtrim(starname,2)               ;trim blanks from starname
    ;online log sheet file
       n_logs=n_elements(logfile)                 ;number of logfiles
       obd0='/mir1/iodspec/'                      ;observation disk
       obd1='/mir1/iodspec/'                 ;observation disk
       obd7='/data/ctio/iodspec/'                 ;observation disk
       obd5='/mir3/iodspec/'                      ;observation disk
endif
       n_logs=n_elements(logfile)                 ;number of logfiles
       if n_logs lt 1 then logfile='/data/ctio/bary/qbcvel.ascii'
       logfile=[logfile]
       bnum=0                                     ;counter
;internal bookeeping struct
        dum={cfstr,obnm:'?',iodnm:'?',bc:0.,z:0.,jd:double(0.),dewar:0,gain:2.5, $
	 cts:long(0),mnvel:0.,mdvel:0.,med_all:0.,errvel:0.,mdchi:0.,nchunk:0, $
	 mdpar:fltarr(20),mnpar:fltarr(20),sp1:0.,sp2:0.,spst:'?',phase:0., $
	 psfpix:[m7_psfpix[0:10], 0.00, 0.00, 0.00,0.00],$
	 psfsig:[m7_psfsig[0:10], 0.00, 0.00, 0.00,0.00]}
	cf=replicate(dum,20000)                     ;  "         "         "
	dum='?'                                    ;string dummy
        obtype='o'
        if starname eq 'iod' then starname='iodine'
        if starname eq 'iodine' then obtype='i'
;This section of code reads the log sheet, creates a list of
;   observations to be analyzed, and makes sure that each
;   observation can be found on disk

for n=0,(n_logs-1) do begin
spawn,'\rm deleteme'
spawnst='grep '+starname+' '+strtrim(logfile(n),2)+ ' >  deleteme'
spawn,spawnst
        close,4 & openr,4,'deleteme'               ;open bary log sheet file
	qq=0                                       ;current log file has not been read
        dum='?'

	while (eof(4) eq 0) do begin               ;begin reading log
           readf,4,dum
;target star?

; obtnum is word location of obtype in logfile, 
;    old style obtnum = 7, new style obtnum = 5
           if strlowcase(strtrim(getwrd(dum,1),2)) eq strlowcase(starname) $
	     and qq eq 0 then begin
	      obtnum=nwrds(dum)-1
	      qq=1
           endif

           if qq gt 0 then if strlowcase(strtrim(getwrd(dum,1),2)) eq strlowcase(starname) and $ 
                strtrim(getwrd(dum,obtnum),2) eq obtype then begin ;if so, begin
               if not keyword_set(noprint) then print,dum                     
               cf(bnum).obnm=strtrim(getwrd(dum,0),2)       ;observation

               if (keyword_set(obnm_filter)) then begin
                   filter_ind = where(obnm_filter eq cf(bnum).obnm)
                   if (filter_ind eq -1) then continue ; not in filter, so skip it.
               endif

	       cf(bnum).bc=float(strtrim(getwrd(dum,2),2))  ;bary correction
	       cf(bnum).jd=double(strtrim(getwrd(dum,3),2)) ;julian date
	       if cf(bnum).jd lt 2440000 then cf(bnum).jd = cf(bnum).jd + 2440000d0
	       cf(bnum).dewar=chip(cf(bnum).obnm,gain)      ;which CCD?
               if cf(bnum).dewar eq -1 then begin           ;CCD not found!
	          print,'Unable to find chip # for observation: ' $
		    +cf(bnum).obnm
                  print,'Assuming Dewar #6 is appropriate.'
		  print,'   Use Control C if you wish to bail out now!'
		  cf(bnum).dewar=6
               endif
	       cf(bnum).gain=gain

               if (keyword_set(dewar)) then begin
                   if (cf(bnum).dewar ne dewar) then continue; not the correct dewar, so skip it.
               endif

	       if strmid(cf(bnum).obnm,0,2) eq 'rh' then obdsk=obd0
	       if strmid(cf(bnum).obnm,0,2) eq 'ra' then obdsk=obd0
	       if strmid(cf(bnum).obnm,0,2) eq 'rb' then obdsk=obd0
	       if strmid(cf(bnum).obnm,0,2) eq 'rd' then obdsk=obd0
	       if strmid(cf(bnum).obnm,0,2) eq 're' then obdsk=obd0
	       if strmid(cf(bnum).obnm,0,2) eq 'rc' then obdsk=obd0
	       if strmid(cf(bnum).obnm,0,2) eq 'rq' then obdsk='/data/ctio/iodspec/'
;	       if strmid(cf(bnum).obnm,0,2) eq 'rs' then obdsk=obd2 df 6-23-99
	       if strmid(cf(bnum).obnm,0,2) eq 'rs' then obdsk='/mir1/iodspec/'
	       if strmid(cf(bnum).obnm,0,2) eq 'rf' then obdsk='/mir1/iodspec/'
	       if strmid(cf(bnum).obnm,0,2) eq 'rg' then obdsk='/mir1/iodspec/'
	       if strmid(cf(bnum).obnm,0,2) eq 'rv' then obdsk=obd4
	       if strmid(cf(bnum).obnm,0,2) eq 'rk' or $
                 strmid(cf(bnum).obnm,0,2) eq 'em'  then begin
                   obdsk = '/mir1/iodspec/'; dave - changed from mir3 so I'm sure to look at my files versus JJ's - 12/2/05.
	       endif
               if strmid(cf(bnum).obnm,0,2) eq 'rj'  then begin
                   obdsk='/mir3/iodspec/'
               endif
	       if strmid(cf(bnum).obnm,0,2) eq 'ru' then begin
		  obdsk=obd3
	       endif
	       if strmid(cf(bnum).obnm,0,2) eq 'rx' then obdsk=obd0
;                           "rx" is post-fix November 1994 Hamilton day-sky
	       dsknm=obdsk+cf(bnum).obnm                    ;Obs. disk name
; The "rk1" observations are interleaved sets of 4 observations, clean up dsknm
	       if strmid(cf(bnum).obnm,0,4) eq 'rk1.' then $
	           dsknm=strmid(dsknm,0,strlen(dsknm)-2)
	       dum=first_el(findfile(dsknm))                ;Obs. on disk?
	       dork=0  &  if dum ne dsknm then dork=1
	       if dork eq 0 then bnum=bnum+1 else begin
;                  if not keyword_set(noprint) then begin
                     if obtype eq 'i' then talk='Iodine: '
                     if obtype eq 'o' then talk='Observation: '
		     talk=talk+dsknm+'  was not found on disk!'
		     print,talk
		     if keyword_set(svlist) then begin
			openu,1,'doplog',/append
			printf,1,talk
			close,1
		     endif
		     print,'   Use Control C if you wish to bail out now!'
;                  endif
               endelse   ;if dork eq 0 then ... else begin
           endif
        end   ;while
        close,4                                          ;Close log file
endfor
        if bnum le 0 then cf=cf(0) else cf=cf(0:bnum-1)  ;Bookeeping structure
;End of the bookeeping/log reading section 

return
end


