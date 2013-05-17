pro rdfts,w,s,wmin,wmax,dfn=dfn,dfd=dfd
;
;  Extract a portion of the FTS Iodine Spectrum
;
;  This routine formerly read /herschel/paul/atlas/ftsiod.bin
;  This routine currently reads /mir/gmarcy/geoff/ftsiod.bin
;
;  INPUTS
;     WMIN   minimum wavelength requested
;     WMAX   maximum wavelength requested
;
;  OUTPUTS
;     W    returned vacuum wavelength scale
;     S    returned FTS I2 spectrum
;
;  KEYWORDS
;     DFN  disk_file_name of fts atlas
;     DFD  disk_file_directory of fts atlas
;            it is assumed that dfd = '/mir1/atlas/' or '/d4/atlas/'
;            unless otherwise specified
;  
;
; Create: Paul Butler/Geoff Marcy, Jan. 1992
; Modified and updated to handle new sacred fts iodine, Oct 10, 1995  PB
; Modified and updated to UVES and other I2 atlases, 6 July 2002 PB
;
;  To convert to "air" wavelengths: use the following line with vactoair
;     wrange=[wmin+1.,wmax+2.]   
;  Find the FTS atlas
;
     if n_params() lt 4 then begin
       print,'Syntax:  rdfts,w,s,wmin,wmax,dfn=dfn,dfd=dfd'
       return
     end

     berkfts='/mir1/atlas/ftsiod.bin'                          ;Berkeley FTS Atlas
     dtmfts ='/mir1/paul/atlas/ftseso50.bin'                   ;DTM FTS Atlas

     if n_elements(dfd) ne 1 then begin          ;FTS disk directory
       dummy=first_el(findfile(berkfts))
       if dummy eq berkfts then dfd='/mir1/atlas/' else begin
         dummy=first_el(findfile(dtmfts))
         if dummy eq dtmfts then dfd='/mir1/paul/atlas/' else begin
            print,'There is no FTS Atlas Directory, such as:'
            print,'   /mir1/atlas/'
            print,'   /d4/atlas/'
            print,'   /mir1/paul/atlas/'
;           print,'   /disks/denali/scratch/paulb/atlas/'
            stop,'A valid FTS Atlas Directory is needed!  Try keyword "dfd".
         endelse
       endelse
     endif

     if n_elements(dfn) ne 1 then dfn='ftslick1.bin'   ;Is the FTS atlas specified?
     ftsfile=strtrim(dfd,2)+strtrim(dfn,2)         ;FTS file
     dummy=first_el(findfile(ftsfile))
     if dummy ne ftsfile then begin                    ;Is "ftsfile" on disk
        print,'Can not find FTS Atlas: '+ftsfile
        stop,'A valid FTS atlas is needed!'
     endif

;ftskeck50_1.35.bin   w0=4989.49557235  n_elements(219553)
;keckfts.pro

;Wavelength zero-point and dispersion for the 1993 FTS Atlas run
;   ftskeck, ftseso, fts1-3a, ftslick
     wav0=20202.018276018d0 
     disp=1.409413936d-2
     maxpix=324746

     case strtrim(dfn,2) of
       'ftslick1.bin': begin          ;New sacred FTS Atlas
           wav0=20202.030866371d0     ;wavenumber of center of pixel 0 in record 0
           disp=1.312212998d-2        ;dispersion
       end

       'ftsiod.bin': begin            ;Old profane FTS Atlas
           wav0=20202.0323d0          ;wavenumber of center of pixel 0 in record 0
           disp=1.40941396d-2         ;extra digit yields better agreement
       end

       'ftsuves_lamp.bin': begin         ;UVES FTS Atlas
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsuves66.bin': begin         ;UVES FTS Atlas
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsuves68.bin': begin         ;UVES FTS Atlas
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsuves70.bin': begin         ;UVES FTS Atlas
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsuves72.bin': begin         ;UVES FTS Atlas
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsuves74.bin': begin         ;UVES FTS Atlas
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas1_lamp.bin': begin       ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.0087737372d0    ;wavenumber of center of pixel 0 in record 0
           disp=1.334083214d-2        ;extra digit yields better agreement
       end

       'ftsas1.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.0087737372d0    ;wavenumber of center of pixel 0 in record 0
           disp=1.334083214d-2        ;extra digit yields better agreement
       end

       'ftsas2.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.0087737372d0    ;wavenumber of center of pixel 0 in record 0
           disp=1.334083214d-2        ;extra digit yields better agreement
       end

       'ftsas3.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.0087737372d0    ;wavenumber of center of pixel 0 in record 0
           disp=1.334083214d-2        ;extra digit yields better agreement
       end

       'ftsas4_lamp.bin': begin       ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
       end

       'ftsas4.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas5.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas6.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas7_lamp.bin': begin       ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas7.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas8.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas9.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftsas10.bin': begin            ;Lick Standard, 10 New Cells, January 2001
           wav0=20202.015444153d0     ;wavenumber of center of pixel 0 in record 0
           disp=6.67041607d-3         ;extra digit yields better agreement
           maxpix=749579   
       end

       'ftslick05.bin': begin         ;Lick Cell, February 2005 
           wav0=20983.301270d0        ;wavenumber of center of pixel 0 in record 0
           disp=0.01334083214d0       ;extra digit yields better agreement
           maxpix=524288   
       end

       'ftsapf05.bin': begin          ;APF Cell same as ftsas??, February 2005 
           wav0=20983.301270d0        ;wavenumber of center of pixel 0 in record 0
           disp=0.01334083214d0       ;extra digit yields better agreement
           maxpix=524288   
       end

       else: begin                    ;All the 1993 FTS Atlas
	   wav0=20202.018276018d0
	   disp=1.409413936d-2
           maxpix=324746
       end
     endcase

     wrange=[wmin-0.6,wmax+0.6]   
     pix=(-1.d8/wrange+wav0)/disp
     pix=[long(pix(0))-1L,long(pix(1))+1L]
     if pix(0) lt 0 then pix(0)=long(0)
     if pix(1) gt long(maxpix) then pix(1)=long(maxpix)
     npix=pix(1)-pix(0)+1L
     w=[1.d8/((dindgen(npix)+pix(0))*(-1.*disp)+wav0)]

     openr,4,ftsfile,/swap_if_little_endian
     a=assoc(4,intarr(npix),2L*pix(0))
     s=double(a(0))
     close,4
     return
end
