pro ck_psf, cdbfiles=cdbfiles

; fischer jan 2013
 num=n_elements(cdbfiles) 
 if num eq 1 then $
 	dop_diagnostics, cdbfile=cdbfiles, /ctio, /rms_map
 	
 if num gt 1 then $
 	dop_diagnostics, cdbfile=cdbfiles[0], /ctio, psfstack=cdbfiles
 	
 	
 stop
 
 end  ;pro		