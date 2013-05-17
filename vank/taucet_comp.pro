pro taucet_comp,hdcopy=hdcopy

restore,'/tous/mir3/vel/vst10700.dat'
cfk=cf3 

restore,'/home/fischer/dop2/vank/vst10700_nov2011.dat'
cfo=cf3

restore,'/home/fischer/dop2/vank/vst10700_jul2012.dat'
cfn=cf3

path='/tous/mir3/iodspec/'

num=n_elements(cfk.obnm)
decker=strarr(num)

for i=0,num-1 do begin
   rdsk,hd,path+cfk[i].obnm,2
   decker[i]=sxpar(hd,'DECKNAME')
endfor

print,decker

if keyword_set(hdcopy) then ps_open,'taucet_comp',/color
plot,sqrt(cfk.cts),cfk.errvel,ps=8,xtitle='!6SNR',ytitle='!6 errvel [m s!u-1!n]',xra=[50,500],/xlog,/xsty,yra=[0,3],/nodata, /ysty
oplot, sqrt(cfk.cts),cfk.errvel, ps=8, col=210
yy=where(strcompress(decker,/rem) eq 'B1')
oplot, sqrt(cfk[yy].cts),cfk[yy].errvel, ps=8, col=155
;y=where(strcompress(decker,/rem) eq 'E3')
;oplot, sqrt(cfk[y].cts),cfk[y].errvel, ps=8, col=190
x=where(cfk.dewar eq 102)
oplot, sqrt(cfk[x].cts),cfk[x].errvel, ps=8, col=90
oplot, sqrt(cfo.cts),cfo.errvel, ps=8, col=222 
;oplot, sqrt(cfn.cts),cfn.errvel, ps=8, col=0
;xyouts,70, 0.6,/data,'!6 CHIRON - post-upgrade',size=1.5
xyouts,100, 2.0,/data,'!6 CHIRON - pre-upgrade',size=1.5
xyouts,260, 2.0,/data,'!6 Keck - pre-upgrade',size=1.5  
xyouts,250, 0.7,/data,'!6 Keck - post-upgrade',size=1.5
plots,[50,500],[1.,1.],linesty=2,col=222
if keyword_set(hdcopy) then ps_close

end 
