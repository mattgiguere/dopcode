pro med_offset,plot=plot

;calculate the median errvels on each side of the CCD change
; determine a weighted (by nobs) median offset
; option to send the patch to the vsts - the old 
; cf3.mnvels will be replaced by cf1.mnvels and then 
; the offset will be applied. 
; fischer june1, 2006

files=findfile('/home/fischer/vank/vstbank/vst*.dat',count=count)

loadct,39
!p.background=255
!p.color=1
           
dif6_13=fltarr(count)  &  n6_13=intarr(count)
dif8_6=fltarr(count)  &  n8_6=intarr(count)
sig6_13=fltarr(count)  & sig8_6=fltarr(count)
openw,1,'offset.txt' 

for i=0,count-1 do begin
    restore,files(i)
    star=strmid(files(i),30,strlen(files(i))-34)
    if keyword_set(plot) then velplot,cf1,/dwr_color
    if stdev(cf1.mnvel) lt 2.5*median(cf1.errvel) then begin
      x1=where(cf1.dewar eq 39,nx1)
      x2=where(cf1.dewar eq 18,nx2)
      x3=where(cf1.dewar eq 24,nx3)
      if nx1 gt 5 and nx2 gt 5 and nx3 gt 5 then begin
        med_d13=median(cf1(x1).mnvel)
        med_d6=median(cf1(x2).mnvel)
        med_d8=median(cf1(x3).mnvel)
        vel6_13 = [cf1(x1).mnvel,cf1(x2).mnvel]
        sig6_13(i)=stdev(vel6_13)
        vel8_6 = [cf1(x2).mnvel,cf1(x3).mnvel]
        sig8_6(i)=stdev(vel8_6)
        dif6_13(i)=med_d6 - med_d13
        n6_13(i)=nx1+nx2

        dif8_6(i)=med_d8 - med_d6
        n8_6(i)=nx2+nx3
;        print,'****************************'
;        printf,1,star, dif6_13(i), dif8_6(i) 
;        print,'Offset 6-13: ',dif6_13(i)
;        print,'Offset 8-6: ',
;        print,'****************************'
        if keyword_set(plot) then begin
          !p.color=1
          plots,[cf1(x1(0)).jd,cf1(x1(nx1-1)).jd],[med_d13,med_d13]
          !p.color=120
          plots,[cf1(x2(0)).jd,cf1(x2(nx2-1)).jd],[med_d6,med_d6]
          !p.color=220
          plots,[cf1(x3(0)).jd,cf1(x3(nx3-1)).jd],[med_d8,med_d8]
        end

      endif
      if nx1 le 5 and nx2 gt 5 and nx3 gt 5 then begin
        med_d6=median(cf1(x2).mnvel)
        med_d8=median(cf1(x3).mnvel)
        vel8_6 = [cf1(x2).mnvel,cf1(x3).mnvel]
        sig8_6(i)=stdev(vel8_6)
        dif8_6(i)=med_d8 - med_d6
        n8_6(i)=nx2+nx3
;        print,'****************************'
;        printf,1,star,dif8_6(i)
;        print,'Offset 8-6: ',dif8_6(i)
;        print,'****************************'
        if keyword_set(plot) then begin
          !p.color=120
          plots,[cf1(x2(0)).jd,cf1(x2(nx2-1)).jd],[med_d6,med_d6]
          !p.color=220
          plots,[cf1(x3(0)).jd,cf1(x3(nx3-1)).jd],[med_d8,med_d8]
        end
      endif
    endif
end

mn6_13=median(dif6_13)
sd6_13=stdev(dif6_13)
mn8_6=median(dif8_6)
sd8_6=stdev(dif8_6)

a1=where(dif6_13 ne 0.0 and dif6_13 gt (mn6_13 - 2.*sd6_13) and dif6_13 lt (mn6_13 +2.*sd6_13) ,na1)
a2=where(dif8_6 ne 0.0 and dif8_6 gt (mn8_6 - 2.*sd8_6) and dif8_6 lt (mn8_6 +2.*sd8_6),na2)


dif6_13=dif6_13(a1)  & n6_13=n6_13(a1)
dif8_6=dif8_6(a2)    & n8_6=n8_6(a2)
sig6_13=sig6_13(a1)  & sig8_6=sig8_6(a2)

ps_open,'offset_med',/color
!p.color=1
plothist,dif8_6,bin=5,/fill,fcolor=220,/fline,forientation=45,$
  xra=[-30,30]
plothist,dif6_13,bin=5,/fill,fcolor=1,/fline,forientation=-45,/overplot
ps_close

err6_13=sig6_13/sqrt(n6_13)
err8_6=sig8_6/sqrt(n8_6)

wt6_13=1/(err6_13)^2.
tot6_13=total(wt6_13)
wt6_13=wt6_13/tot6_13

wt8_6=1/(err8_6)^2.
tot8_6=total(wt8_6)
wt8_6=wt8_6/tot8_6

wm_6_13 = total(dif6_13*wt6_13)
wm_8_6 = total(dif8_6*wt8_6)

form1='(a10,i6,f9.2)'
openw,2,'woffset.txt'
printf,2,'Dewar 6 to 8 Transition '
for i=0,na1-1 do printf,2,strmid(files(a1(i)),30,strlen(files(a1(i)))-34),n6_13(i),dif6_13(i),format=form1
printf,2,'Dewar 6 to 8 Transition '
for i=0,na2-1 do printf,2,strmid(files(a2(i)),30,strlen(files(a2(i)))-34),n8_6(i),dif8_6(i),format=form1
close,2

ps_open,'hist_offset'
plothist,dif8_6,bin=2.,yra=[0,12],/ysty,xtit='!6 RV Offset from Dewar 6 to Dewar 8 [m s!u-1!n]'
ps_close

stop
print,'Weighted mean Dewar 6-13: ',wm_6_13
print,'Weighted mean Dewar 8-6: ',wm_8_6
wm_6_13 = -1.0*wm_6_13
wm_8_6 = -1.0*wm_8_6

ans=''
read,'Make a patch to the vsts? (y/n): ',ans
if ans eq 'y' then patch_cf3,wm_6_13,wm_8_6,/redo

close,1
stop


end

