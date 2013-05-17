pro ck_outlier, cdb=cdb,fix=fix

if ~keyword_set(cdb) then cdb='cdbbHR4743_rqa05.3054'
path='/mir7/files/'

restore,path+cdb

; check continuity of the wavelength soln
; j = order
; i = chunk in that order

ord=chunk_arr.ordt
i=uniq(ord)
ord=ord[i] 
nord=n_elements(ord) 

!p.multi=[0,2,1]
for j=0, nord-1 do begin 
  x=where(chunk_arr.ordt eq ord[j],nchnks)
  wav=dblarr(nchnks) 
  disp=dblarr(nchnks) 
  chnks=chunk_arr[x]   ; just that order 
  for i=0, nchnks-1 do wav[i]=(*chnks[i].free_par)[2].wcof[0]
  for i=0, nchnks-1 do disp[i]=(*chnks[i].free_par)[2].wcof[1]
  plot,wav,/xsty,/ysty,charsize=1.7,ps=8
  plot,disp,/xsty,/ysty,charsize=1.7,ps=8
  coef=poly_fit(findgen(nchnks),disp,1)
  ndisp=poly(findgen(nchnks),coef)
  oplot,ndisp,col=222,ps=8
  if keyword_set(fix) then for i=0,nchnks-1 do (*chunk_arr[x[i]].free_par)[2].wcof[1]=ndisp[i]
end

  if keyword_set(fix) then save,chunk_arr,f=path+cdb

end



