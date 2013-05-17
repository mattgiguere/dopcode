pro conv_vd2cb, dopenv, tag=tag, vdiodnm=vdiodnm, obsnm=obsnm, observatory=observatory


; after bootstrapping into a vd with dop_bearing_het.pro
; this program will format into the structures 
; needed by dop_main.
; fischer, 8/4/08  sfsu

  nchks=dopenv.n_chunks
  path=dopenv.files_dir

  if ~keyword_set(obsnm) then read,'Enter obsnm: ',obsnm

  if ~keyword_set(tag) then tag='c'
  outfile='cd'+tag+'b'+dopenv.obj_nm+'_'+obsnm
  
  restore,dopenv.files_dir+vdiodnm          ;vd

  chunk_arr=dop_chunk_replicate(nchks, dopenv)
  nchunks=n_elements(vd.ordt)
  chunk_arr.ordt = vd.ordt
  chunk_arr.ordob = vd.ordt
  chunk_arr.pixt = vd.pixt
  chunk_arr.pixob = vd.pixt
  chunk_arr.weight = 1.0
  init_z = 0
 
  for i=0,nchunks-1 do begin
     (*chunk_arr[i].free_par)[*].wcof[0] = vd[i].wcof[0] + vd[i].w0
     (*chunk_arr[i].free_par)[*].wcof[1] = vd[i].wcof[1]
     (*chunk_arr[i].free_par)[*].amp[0] = dopenv.psfsig[0]
     for j=1,10 do (*chunk_arr[i].free_par)[0].amp[j] = vd[0].par[j]
     (*chunk_arr[i].free_par)[*].z=init_z
  endfor

     save,chunk_arr,f=path+outfile
stop
  end
