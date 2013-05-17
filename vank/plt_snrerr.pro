pro plt_snrerr, hdcopy=hdcopy

;fischer jan29,2013

; single observations, snr~160 
vank,'10700','adg','vdadg',cf1,cf3,/ctio,dfd='/tous/mir7/files_df/t1/'
restore,'vstbank/vst10700.dat'  
velplot,cf3,4./24.,d,v,e,n,c,dum1
print,'SNR: ',sqrt(dum1.cts) 
print,'Err: ',mean(dum1.errvel)
print,'RMS: ',stddev(dum1.mnvel) 

;double observations, snr ~224
vank,'10700','adg','vdadg',cf1,cf3,/ctio,dfd='/tous/mir7/files_df/t2/'
restore,'vstbank/vst10700.dat'  
velplot,cf3,4./24.,d,v,e,n,c,dum2
print,'SNR: ',sqrt(dum2.cts) 
print,'Err: ',mean(dum2.errvel)
print,'RMS: ',stddev(dum2.mnvel) 


;triple observations, snr ~275
vank,'10700','adg','vdadg',cf1,cf3,/ctio,dfd='/tous/mir7/files_df/t3/'
restore,'vstbank/vst10700.dat'  
velplot,cf3,4./24.,d,v,e,n,c,dum3
print,'SNR: ',sqrt(dum3.cts) 
print,'Err: ',mean(dum3.errvel)
print,'RMS: ',stddev(dum3.mnvel) 


; all observations, snr~400 
t0=16112.
vank,'10700','adg','vdadg',cf1,cf3,/ctio,dfd='/tous/mir7/files_df/'
restore,'vstbank/vst10700.dat'  
x=where(cf3.jd gt 16115. and cf3.jd lt 16127.,nx) 
cf=cf3[x]
velplot,cf,4./24.,d,v,e,n,c,dumall
print,'SNR: ',sqrt(dumall.cts) 
print,'Err: ',mean(dumall.errvel)
print,'RMS: ',stddev(dumall.mnvel) 

; data plots 
!x.margin=[15,10]
!y.margin=[8,6]
if keyword_set(hdcopy) then !p.font=1 else !p.font=-1
!p.charsize=1.5
!p.color=0

if keyword_set(hdcopy) then ps_open,'taucet_2wks',/color
plot,dumall.jd-t0,dumall.mnvel,ps=8,symsize=0.4, yra=[-20,20],/ysty,xra=[3.,16.],/xsty, $
     xtitle='!6 JD-2456110',ytitle='!6 Radial Velocity [m s!u-1!d]'
oploterr,dumall.jd-t0,dumall.mnvel,dumall.errvel,1
str0 = '!6Tau Ceti, Orders 13 - 33 '
xyouts,10, -12,/data,str0
xyouts,10, -14,/data,'!6 CTIO / CHIRON data'
str1 = Greek('sigma')+'!6!dINT!n = '+strmid(strtrim(mean(dumall.errvel),2),0,4)+' ms!u-1!n'
str2='!6 RMS ='+strmid(strtrim(stddev(dumall.mnvel),2),0,4)+' ms!u-1!n'
xyouts, 11,17,/data,str2,size=1.5
xyouts, 11,15,/data,str1,size=1.5
if keyword_set(hdcopy) then ps_close

; impact of fewer orders
!p.color=0
t0=16112.
ordt=indgen(15)+13
vank,'10700','adg','vdadg',cf1,cf3,/ctio,dfd='/tous/mir7/files_df/',order=ordt
restore,'vstbank/vst10700_ordr.dat'  
x=where(cf3.jd gt 16115. and cf3.jd lt 16127.,nx) 
cf=cf3[x]
velplot,cf,4./24.,d,v,e,n,c,dumord
if keyword_set(hdcopy) then !p.font=1 else !p.font=-1
!p.charsize=1.5
!p.color=0
if keyword_set(hdcopy) then ps_open,'taucet_2wks_red_ord',/color
plot,dumord.jd-t0,dumord.mnvel,ps=8,symsize=0.4,yra=[-20,20],xra=[3,16],/xsty,/ysty
oploterr,dumord.jd-t0,dumord.mnvel,dumord.errvel,1   ;psym=3
str0 = '!6Tau Ceti, Orders 13 - 27'
xyouts,10, -12,/data,str0
xyouts,10, -14,/data,'!6 CTIO / CHIRON data'
str1 = Greek('sigma')+'!6!dINT!n = '+strmid(strtrim(mean(dumord.errvel),2),0,4)+' ms!u-1!n'
str2='!6 RMS ='+strmid(strtrim(stddev(dumord.mnvel),2),0,4)+' ms!u-1!n'
xyouts, 11,17,/data,str2,size=1.5
xyouts, 11,15,/data,str1,size=1.5
if keyword_set(hdcopy) then ps_close


; how do the errror and the RMS scale with SNR? 
snr_arr=[sqrt(mean(dum1.cts)), sqrt(mean(dum2.cts)), sqrt(mean(dum3.cts)), sqrt(mean(dumall.cts))]
err_arr=[mean(dum1.errvel), mean(dum2.errvel), mean(dum3.errvel), mean(dumall.errvel)]
rms_arr=[stddev(dum1.mnvel), stddev(dum2.mnvel), stddev(dum3.mnvel), stddev(dumall.mnvel)]
coef1=poly_fit(snr_arr, err_arr, 2)
coef2=poly_fit(snr_arr, rms_arr, 2)
xarr=findgen(24)*10. + 160.

if keyword_set(hdcopy) then ps_open,'snr_err',/color
plot,snr_arr,err_arr,ps=8,yra=[0,2],xtitl='!6 Signal-to-noise ratio',xra=[150,400],/xsty,ytitl='!6 [m s!u-1!n]' 
oplot,snr_arr,rms_arr,ps=8,col=222  
oplot,xarr,poly(xarr,coef1),linest=2  
oplot,xarr,poly(xarr,coef2),linest=2,col=222  
xyouts,300, 0.2,/data,'!6CTIO / CHIRON Data' 
!p.color=222
xyouts,270, 1.62,/data,'!6RMS scatter'     
!p.color=0
xyouts,270, 1.5,/data,'!6Single measurement errors',size=1.5   
if keyword_set(hdcopy) then ps_close

end ;pro
