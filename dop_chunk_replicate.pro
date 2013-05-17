FUNCTION DOP_CHUNK_REPLICATE, nchunk, dopenv

  npar_psf=n_elements(dopenv.psfpix)
  for i=0,nchunk-1 do begin 

 ; DEFINE THE CHUNK STRUCTURE TAG NAMES
 ; WITH INFORMATION SPECIFIC TO EACH CHUNK
      one_chunk = {ordt:0,                    $ ;order of template  
            ordob:0,                    $ ;order of obs
            pixt:0,                     $ ;starting chunk pixel, template
            pixob:0,                    $ ;starting chunk pixel, obs
            cts:long(0),                $ ;photon counts 
            sobs:ptr_new(fltarr(dopenv.n_pix)),$ ;pointer to observed, FTS and template spec
            siod:ptr_new(/allocate),    $ ;pointer for iodine spectrum
            wiod:ptr_new(/allocate),    $ ;pointer for iod wavelengths
            sis:ptr_new(/allocate),     $ ;flux of the intrinsic spectum
            wis:ptr_new(/allocate),     $ ;wavelength of the intrinsic spectrum
            smod:ptr_new(/allocate),    $ ;model spectrum (wavelengths in free_par)
            free_par: ptr_new(replicate({wcof:dblarr(dopenv.n_wcof), z:0.0d, $
                                         offset:1.0d, amp:dblarr(npar_psf)},dopenv.n_iter)), $
            free_ind:intarr(npar_psf+4,dopenv.n_iter), $ ;1 if free, 0 if not
            fit:fltarr(dopenv.n_iter),    $
            ip:fltarr(121,dopenv.n_iter), $ ;store IP
            weight:0.0,                   $ ;chunk weight (from ISS) 
            gdpix: 0}                     ;number of good pixels (filter known bad pixels)

      if i eq 0 then chunk_arr=one_chunk else chunk_arr=[chunk_arr, one_chunk]
   endfor

   return, chunk_arr

end
