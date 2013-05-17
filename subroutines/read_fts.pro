pro read_fts, w1, w2, dopenv, wiod, siod

; fischer Jan27, 2013
; INPUT
; w1, w2 (wavelength range boundaries)
; dopenv structure with basic info, including paths

; OUTPUT
; wiod, siod


; different calls for different atlases...

; Keck or Lick 
  if dopenv.fts_atlas eq 'ftsas4.bin' or dopenv.fts_atlas eq 'ftskeck50.bin' or $
     dopenv.fts_atlas eq 'ftseso50.bin' or dopenv.fts_atlas eq 'ftslick05.bin' then $
     rdfts,wiod,siod,w1,w2,dfn=dopenv.fts_atlas, dfd=dopenv.fts_dir

; (too) highly sampled CHIRON
;    if dopenv.fts_atlas eq 'iodine_ctiov1f_pnnl.dat' then $
;       read_iodine,'ctiov1f','pnnl',40.0, wav[pix1]-ipad,wav[pix2]+ipad,wiod,siod,nadd=1 

; CHIRON
  if dopenv.fts_atlas eq 'ctiov1f_pnnl_red5.dat' then begin
     read_iodine,'ctiov1f','pnnl',40.0,w1,w2,wiod,siod,nadd=0
  endif 
            
  if dopenv.fts_atlas eq 'iodine_keck_pnnl.sav' then $
     read_iodine,'keck','pnnl',63.0,w1,w2,wiod,siod,nadd=1 
            
  if dopenv.fts_atlas eq 'ftshet.dat' then $ 
     rdi2,wiod,siod,w1, w2,dopenv
            
  if dopenv.fts_atlas eq 'iodine_p4_pnnl_lo_apod.dat' then $
     read_iod_p4,dopenv.fts_atlas,w,s,w1,w2,pad=1.

end



