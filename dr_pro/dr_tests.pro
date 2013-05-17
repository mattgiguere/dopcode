pro dr_tests, test1=test1, test2=test2, run=run

; fischer may 17, 2010 
;  checks of the alpha cen data


; the night of May 8,9,10, several sets of iodine observations were
; obtained to check PSF stability through the night 

if keyword_set(test1) then begin
  ; TEST 1: PSF stability within one night and from night to night for
  ; the iodine observations for May 8, 9, 10
    restore,'/mir7/bary/qbcvel.dat'
    xn1=where(bcat.objnm eq 'iodine' and strmid(bcat.obsnm,0,5) eq run and $
              fix(strmid(bcat.obsnm,6,4)) gt 5150 and fix(strmid(bcat.obsnm,6,4)) lt 5660,nxn1)
    xn2=where(bcat.objnm eq 'iodine' and strmid(bcat.obsnm,0,5) eq run and $
              fix(strmid(bcat.obsnm,6,4)) gt 5700 and fix(strmid(bcat.obsnm,6,4)) lt 6310,nxn2)
    xn3=where(bcat.objnm eq 'iodine' and strmid(bcat.obsnm,0,5) eq run and $
              fix(strmid(bcat.obsnm,6,4)) gt 6400 and fix(strmid(bcat.obsnm,6,4)) lt 6810,nxn3)

    ind=[xn1,xn2,xn3]
    nstars=nxn1+nxn2+nxn3

    restore,'/mir7/files/vdbiod_'+bcat[ind[0]].obsnm
    nip=n_elements(vd[0].psf)
    nchunks=n_elements(vd.ordt)
    psfarr=fltarr(nip,nstars,nchunks)
    chisq=fltarr(nstars) 
    badflag=intarr(nstars) 

    for i=0,nstars-1 do begin
       vdnm='vdbiod_'+bcat[ind[i]].obsnm
       restore,'/mir7/files/'+vdnm
       chisq[i]=median(vd.fit)
       for j=0,nchunks-1 do begin
          psfarr[*,i,j]=vd[j].psf[*]
       endfor
    endfor

    medpsf=fltarr(nip,nchunks)  ; median psf for all observations
    sigpsf=fltarr(nip,nchunks)  ; sig psf for all obs
  for j=0,nchunks-1 do begin
     for l=0,nip-1 do medpsf[l,j]=median(psfarr[l,*,j])
     for l=0,nip-1 do sigpsf[l,j]=stddev(psfarr[l,*,j])
  endfor

  for j=0,nchunks-1 do begin     
     plot,medpsf[*,j],col=1,/xsty, thick=3,yra=[-0.05,0.7]
     for i=0,nstars - 1 do begin
        sind=where(abs(medpsf[*,j] - psfarr[*,i,j]) gt 4*sigpsf[*,j] and abs(medpsf[*,j] - psfarr[*,i,j]) gt 0.02,nsind)
        if nsind gt 1 then badflag[i]=badflag[i]+1 
        if nsind gt 1 then print,bcat[ind[i]],nsind,medpsf[sind,j] - psfarr[sind,i,j]
        if nsind gt 1 then oplot,psfarr[*,i,j],col=222 $
        else oplot,psfarr[*,i,j],col=i
     endfor
     oplot,medpsf[*,j],col=1,thick=3
     str1='!6 RMS: '+strmid(strcompress(string(mean(sigpsf[where(medpsf[*,j] ne 0),j])),/rem),0,5)
     xyouts,80,0.7*max(medpsf[*,j]),/data,str1
  endfor

  x=where(badflag gt 1,nx)
  if nx gt 0 then print,ff[x], badflag[x]

endif  ;TEST 1

if keyword_set(test2) then begin
  ; TEST 2: PSF stability within one night and from night to night for
  ; stellar observations on May 8, 9, 10
    restore,'/mir7/bary/qbcvel.dat'
    xn1=where(bcat.objnm eq '128620' and strmid(bcat.obsnm,0,5) eq run and $
              fix(strmid(bcat.obsnm,6,4)) gt 5200 and fix(strmid(bcat.obsnm,6,4)) lt 5540,nxn1) ; 5660,nxn1)
    xn2=where(bcat.objnm eq '128620' and strmid(bcat.obsnm,0,5) eq run and $
              fix(strmid(bcat.obsnm,6,4)) gt 5700 and fix(strmid(bcat.obsnm,6,4)) lt 6310,nxn2)
    xn3=where(bcat.objnm eq '128620' and strmid(bcat.obsnm,0,5) eq run and $
              fix(strmid(bcat.obsnm,6,4)) gt 6400 and fix(strmid(bcat.obsnm,6,4)) lt 6880,nxn3)

    ind=[xn2] ;,xn2]
    nstars=nxn2 ;+nxn2

    restore,'/mir7/files/vdb128620_'+bcat[ind[0]].obsnm
    nip=n_elements(vd[0].psf)
    nchunks=n_elements(vd.ordt)
    psfarr=fltarr(nip,nstars,nchunks)
    chisq=fltarr(nstars) 
    badflag=intarr(nstars) 

    for i=0,nstars-1 do begin
       vdnm='vdb128620_'+bcat[ind[i]].obsnm
       restore,'/mir7/files/'+vdnm
       chisq[i]=median(vd.fit)
       for j=0,nchunks-1 do begin
          psfarr[*,i,j]=vd[j].psf[*]
       endfor
    endfor

    medpsf=fltarr(nip,nchunks)  ; median psf for all observations
    sigpsf=fltarr(nip,nchunks)  ; sig psf for all obs
  for j=0,nchunks-1 do begin
     for l=0,nip-1 do medpsf[l,j]=median(psfarr[l,*,j])
     for l=0,nip-1 do sigpsf[l,j]=stddev(psfarr[l,*,j])
  endfor

  for j=0,nchunks-1 do begin
     plot,medpsf[*,j],col=1,/xsty, thick=3,yra=[-0.05,0.7]
     for i=0,nstars - 1 do begin
        sind=where(abs(medpsf[*,j] - psfarr[*,i,j]) gt 4*sigpsf[*,j] and abs(medpsf[*,j] - psfarr[*,i,j]) gt 0.02,nsind)
        if nsind gt 1 then badflag[i]=badflag[i]+1 
        if nsind gt 1 then print,bcat[ind[i]],nsind,medpsf[sind,j] - psfarr[sind,i,j]
        if nsind gt 1 then oplot,psfarr[*,i,j],col=222 $
        else oplot,psfarr[*,i,j],col=i
     endfor
     oplot,medpsf[*,j],col=1,thick=3
     str1='!6 RMS: '+strmid(strcompress(string(mean(sigpsf[where(medpsf[*,j] ne 0),j])),/rem),0,5)
     xyouts,80,0.7*max(medpsf[*,j]),/data,str1
  endfor

  x=where(badflag gt 1,nx)
  if nx gt 0 then print,ff[x], badflag[x]
stop
 endif ;TEST2
    


stop
end
