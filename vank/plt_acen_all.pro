pro plt_acen_all,star=star


restore,'vst'+star+'.dat'  ; 2009
cfo=cf3

clean_acen,tag='adg',star=star
restore,'vstbank/vst'+star+'_clean_adg.dat'
cfn=cf3
if star eq 128620 then sft=528.
cfn.mnvel=cfn.mnvel+sft

cf=[cfo,cfn]
!x.margin=[15,10]
!y.margin=[6,6]
!p.charsize=1.6
!x.charsize=1.6
!y.charsize=1.6
coef=poly_fit(cfo.jd,cfo.mnvel,1)
plot,cf.jd-14500.,cf.mnvel,ps=8,symsize=0.2,xtitl='!6JD-2454500',$
	ytitle='!6 RV',yra=[-200,1000],/ysty,title='!6 HD 128620'
oplot,cf.jd-14500.,poly(cf.jd,coef),col=222

ans='y'
read,'Adjust the offset for the new data? (y/n) ',ans
if ans eq 'y' then begin
  repeat begin
	shft=0.0
	read,'Enter the shift for the new velocity set [m/s] ',shft
	cfn.mnvel=cfn.mnvel+shft
	cf=[cfo,cfn]
	plot,cf.jd-14500.,cf.mnvel,ps=8,symsize=0.2,xtitl='!6JD-2450000',ytitle='!6 RV',$
		yra=[-200,1000],/ysty,title='!6 HD 128620'
	oplot,cf.jd-14500.,poly(cf.jd,coef),col=222
	read,'Adjust the offset for the new data? (y/n) ',ans
  endrep until ans eq 'n'
endif

xyouts,14620.-14500.,-50,/data,'!6 2008',size=1.7
xyouts,14850.-14500.,50,/data,'!6 2009', size=1.7
xyouts,16080.-14500.,580,/data,'!6 2012', size=1.7

stop
end ; pro 