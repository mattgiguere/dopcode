FUNCTION DOP_CHUNK_UPDATE, chunk, parinfo, ip, syn_fit=syn_fit, chi=chi, pass=pass, iod_flag=iod_flag

   ; FILL THE CHUNKS WITH THE MODEL PARAMETERS AFTER PASS 1
    c_light = 2.99792458d8

if pass eq 1 then begin

    names = ['wcof0', 'wcof1','z','offset','amp0','amp1','amp2', $
          'amp3', 'amp4','amp5','amp6','amp7','amp8', 'amp9', 'amp10', $
          'amp11', 'amp12','amp13','amp14','amp15','amp16'] 
    npar = n_elements(parinfo.value)
    chunk.fit[pass] = chi  
    chunk.ip[*,pass]=ip
    (*chunk.smod)=syn_fit

  ; UPDATE FREE PARAMETERS
    (*chunk.free_par)[pass].wcof[0] = parinfo[0].value
    (*chunk.free_par)[pass].wcof[1] = parinfo[1].value
    (*chunk.free_par)[pass].z = parinfo[2].value/c_light
    (*chunk.free_par)[pass].offset = parinfo[3].value
    (*chunk.free_par)[pass].amp[0] = parinfo[4].value
    (*chunk.free_par)[pass].amp[1] = parinfo[5].value
    (*chunk.free_par)[pass].amp[2] = parinfo[6].value
    (*chunk.free_par)[pass].amp[3] = parinfo[7].value
    (*chunk.free_par)[pass].amp[4] = parinfo[8].value
    (*chunk.free_par)[pass].amp[5] = parinfo[9].value
    (*chunk.free_par)[pass].amp[6] = parinfo[10].value
    (*chunk.free_par)[pass].amp[7] = parinfo[11].value
    (*chunk.free_par)[pass].amp[8] = parinfo[12].value
    (*chunk.free_par)[pass].amp[9] = parinfo[13].value
    (*chunk.free_par)[pass].amp[10] = parinfo[14].value
    (*chunk.free_par)[pass].amp[11] = parinfo[15].value
    (*chunk.free_par)[pass].amp[12] = parinfo[16].value
    (*chunk.free_par)[pass].amp[13] = parinfo[17].value
    (*chunk.free_par)[pass].amp[14] = parinfo[18].value
    (*chunk.free_par)[pass].amp[15] = parinfo[19].value
    (*chunk.free_par)[pass].amp[16] = parinfo[20].value
 endif  ;pass 1 

if pass eq 2 then begin
    names = ['wcof0', 'wcof1','z','offset']
    npar = n_elements(names)

    chunk.fit[pass] = chi  
    chunk.ip[*,pass]=ip
    (*chunk.smod)=syn_fit

  ; UPDATE FREE PARAMETERS
    (*chunk.free_par)[pass].wcof[0] = parinfo[0].value
    (*chunk.free_par)[pass].wcof[1] = parinfo[1].value
    (*chunk.free_par)[pass].z = parinfo[2].value/c_light
    (*chunk.free_par)[pass].offset = parinfo[3].value

  ; PSF TAKEN FROM PASS 1, NO FREE PSF PAR'S FOR PASS 2
    fp=1  ; first Doppler pass par's
    for i=0,16 do (*chunk.free_par)[pass].amp[i]= 0.0   ;(*chunk.free_par)[fp].amp[i]
 endif  ; pass 2

if pass eq 3 then begin
    names = ['wcof0', 'wcof1','z','offset']
    npar = n_elements(names)

    chunk.fit[pass-1] = chi  
    chunk.ip[*,pass-1]=ip
    (*chunk.smod)=syn_fit

  ; UPDATE FREE PARAMETERS
    (*chunk.free_par)[pass-1].wcof[0] = parinfo[0].value
    (*chunk.free_par)[pass-1].wcof[1] = parinfo[1].value
    (*chunk.free_par)[pass-1].z = parinfo[2].value/c_light
    (*chunk.free_par)[pass-1].offset = parinfo[3].value

  ; PSF TAKEN FROM PASS 1, NO FREE PSF PAR'S FOR PASS 2
    fp=1  ; first Doppler pass par's
    for i=0,16 do (*chunk.free_par)[pass-1].amp[i]= 0.0   ;(*chunk.free_par)[fp].amp[i]
 endif  ;pass 3


  mod_chunk=chunk

  return, mod_chunk

end
