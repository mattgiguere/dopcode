pro cf,starname,barydsst,cfout,cfiod,logfile=logfile, tape=tape, $
       first_tape=first_tape,noi2=noi2,obdsk=obdsk,obnm_filter=obnm_filter,$
       dewar=dewar
;This code drives the PSF and velocity analysis
;
;INPUT:
;   starname (input string)   Examples: '509' or '4983' or 'GL897'
;   bardsst  (input scalar)   Barycentric correction of the template
;OUTPUT:
;   cfout    (structure)       Information about each observation
;
;OPTIONAL but highly RECOMMENDED:
;   logfile (keyword string)  Specifies log-sheet file to be searched 
;			      rather than the default (bcvel.ascii)
;                             example:  logfile='/mir1/bary/bcvel.cep'
;                             logfile=['/mir1/bary/bcvel.cep',/mir1/bary/bcvel.ascii']
;OPTIONAL:
;   cfiod    (structure)      Information about "Iodine" observations
;   first_tape (keyword string)
;                             Specifies that information is wanted
;                             only for observations starting with
;                             a particular tape, and all subsequent
;                             observations
;                             example:  first_tape='rb02'
;   tape     (keyword string or string array)
;			      Specifies a tape or set of tapes
;			      example:  tape='ra40'
;                             example:  tape=['rh50','rh51','rh52','rh53']
;   noi2     (keyword - on/off)  if "/noi2" invoked, the code will not
;		              bother getting I2 information associated
;		              with observations
;
;
;Created August 28, 1993  R.P.B.
;Modified May 1994
;Millions of comments added August 1995
IF n_params () lt 3 then begin
  print,'-------------------------------------------------------------------'
  print,' SYNTAX:'
  print,' '
  print,' IDL> cf, starname, barydsst, cfout, cfiod, logfile=logfile, tape=tape, $'
  print,'          first_tape=first_tape,noi2=noi2'
  print,' '
  print,'EXAMPLE:'
  print,'IDL>  starname=''8086'' '
  print,'IDL>  barydsst=1228.3   barycentric correction associated with DSST'
  print,'NB:  It is recommended that keyword, LOGFILE, be specified.'
  print,' logfile (keyword string)  Specifies log-sheet file(s) to be searched'
  print,'                           rather than the default (bcvel.ascii)'
  print,'IDL>  logfile=''/mir1/bary/bcvel.cep'' '
  print,"IDL>  logfile=['/mir1/bary/bcvel.cep','/mir1/bary/bcvel.ascii']"
  print,'cfout (output structure)   Information about each observation'
  print,'cfiod (optional output structure) Information about "Iodine" observations
  print,'tape  (keyword string or string array)'
  print,'                           Specifies a tape or set of tapes
  print,'IDL>  tape=''ra16''       '
  print,"IDL>  tape=['rh50','rh51','rh52','rh53']"

  print,'IDL>  cf,starname,barydsst,cfout,logfile=logfile,tape=tape'
  print,"IDL>  save,file='/mir1/paul/idle/dop/cf509_rb02.dat',cfout"
  print,' '
  print,'-----------------------------------------------------------------'
  RETURN
ENDIF
  ;
       c = 2.99792458d8                           ;speed of who?
       bnum=0                                     ;counter
       starname=strtrim(starname,2)               ;trim blanks from starname
       if n_elements(tape) lt 1 then tape =['?'] else tape=[tape]

       fd0='/mir1/bary/qbcvel.ascii'
       fd7='/data/ctio/bary/qbcvel.ascii'
       fd1='/mir7/bary/bcvel.cep'
       ob0='/data/ctio/iodspec/'
       ob1='/mir1/paul/cepspec/'


;This section of code reads the log sheet, creates a list of
;   observations to be analyzed, and makes sure that each
;   observation can be found on disk

;reading iodines and observations from logsheets
    if n_elements(logfile) lt 1 then logfile=[fd7] 

;    if not keyword_set(noi2) then $
;       rdbcvel,'iodine',cfiod,obdsk,logfile=logfile,/noprint,obdsk=obdsk,obnm_filter=obnm_filter,dewar=dewar
    obdsk='/data/ctio/iodspec/'


    rdbcvel,strlowcase(starname),cf,obdsk,logfile=logfile,/noprint,obdsk=obdsk,obnm_filter=obnm_filter,dewar=dewar

    dum=where(cf.obnm ne '?',bnum)
    cf=cf(dum)
;end reading iodines and observations from logsheets
        if n_elements(cfiod) gt 1 then begin
	   dum=sort(cfiod.jd)  &  cfiod=cfiod(dum)
        endif
        if n_elements(cf) gt 1 then begin
           dum=sort(cf.jd)  &  cf=cf(dum)
        endif
	cf.z=(barydsst - cf.bc)/c                        ;Doppler Z guess

;is the tape keyword called?
	if tape(0) ne '?' then for n=0,bnum-1 do begin
	     tpnm=strmid(cf(n).obnm,0,5) ;tape name
             tcase=0                                     ;initialize tcase 
;If tcase eq 1, tape name matches tape keyword. If tcase eq 0, no match.
	     for q=0,(n_elements(tape)-1) do if tpnm eq tape(q) then tcase=1
	     if tcase eq 0 then cf(n).dewar=-1           ;no tape name match
        endfor ;tape(0) eq '?'                           ;end tape keyword
	dum=where(cf.dewar ne -1,bnum)
	if bnum gt 0 then cf=cf(dum) else begin
	   print,'No matching observations for star: '+starname
	   print,'For tapes: '+tape
	   retall
        endelse

        if n_elements(first_tape) eq 1 then begin
           first_tape=strtrim(first_tape,2)
	   cfdum=cf
           if not keyword_set(noi2) then cfdum=cfiod
           dum=where(strmid(cfdum.obnm,0,5) eq first_tape,ndum)
	   if ndum lt 1 then begin
	      cfdum=cf
              dum=where(strmid(cfdum.obnm,0,5) eq first_tape,ndum)
           endif
           if ndum ge 1 then begin
             juldt=cfdum(dum(0)).jd-.6
             dum=where(cf.jd gt juldt)
             cf=cf(dum)
           endif
           bnum=n_elements(cf)
print,'bnum',bnum
           if bnum lt 1 then begin
              talk='No post- '+first_tape+' observations of '+starname
              talk=talk+' were found on disk!'
              print,talk
              retall
           endif
        endif   ;keyword first_tape

if n_elements(tape) eq 1 then tape=tape(0)

cfout=cf

;End of the bookeeping/log reading section 
return
end
