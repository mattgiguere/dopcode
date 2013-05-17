pro str_128620,catName=catName

; create or update an existing structure 
; code by John Brewer 
; modified for doppler code by D. Fischer 22 Dec 2009

restore,'/mir7/bary/qbcvel.dat'
x=where(bcat.objnm eq '128620' and bcat.obtype eq 'o',nobs)
bcat=bcat[x]
nobs=n_elements(bcat)

  catName='db_128620'
  catPath='/home/fischer/dop2/'
; make a structure for batch_dop 
  catTmpl={name:'128620',dir:'/home/fischer/dop2/dop_results/',runnum:'',obsnm:'',$
       status:'NOT RUN',runstartdate:'',runenddate:'',runnotes:''}
  cat=replicate(catTmpl,nobs) 

 ; check for an existing catalog file and catalog lock
  lockInfo = File_Info(catPath+catName+".lock")
  waitCount = 0
  while (lockInfo.Exists and waitCount le 15) do begin
     print,"Waiting for catalog lock to be released...("+strtrim(string(waitCount*2),2)+")"
     wait, 2 
     lockInfo = File_Info(catPath+catName+".lock")
     waitCount = waitCount + 1
  end 
  if (lockInfo.Exists) then begin
     print,"Unable to obtain lock on catalog file after ",strtrim(string(waitCount*2),2)+"seconds."
     stop
  endif

  ; Lock the catalog 
  openw,catlun, catPath+catName+".lock",/get_lun
  printf,catlun, "Writing new stars to catalog (db_128620.dat)"
  free_lun,catlun

  ; Edit the catalog if it exists 
  catInfo = File_Info(catPath+catName+".dat")
  if (catInfo.Exists) then restore, catPath+catName+".dat"  

  for i=0L,nobs-1 do begin
     x=where(bcat[i].obsnm eq cat.obsnm, nx) 
     if nx eq 1 then cat[i]=cat[x]
     if nx eq 0 then begin      ; not in cat
        cat[i].obsnm = bcat[i].obsnm
        cat[i].runnum=strmid(cat[i].obsnm,0,5)
        ff=findfile('/mir7/files/cdbb128620_'+cat[i].obsnm,count=count)
        if count eq 1 then cat[i].status = 'COMPLETED'
     endif
  endfor
  
  save,cat,f=catPath+catName+".dat"
  
 ; remove the lock file
  File_Delete, catPath+catName+".lock"

end

