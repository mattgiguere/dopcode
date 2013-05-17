pro comp_vst

; compare velocities created with dsst128620b_rqa10.dat (JJ's deconv code)
; with velocities created with dsst128620m_rqa10.dat (solar morph code) 

restore,'vstbank/vst128620m.dat'
cfm=cf3

restore,'vstbank/vst128620b.dat'
cfb=cf3 

num=n_elements(cfb)
chiarrb=fltarr(num)
errarrb=fltarr(num)
chiarrm=fltarr(num)
errarrm=fltarr(num)

for i=0,num-1 do begin
   x=where(cfb[i].obnm eq cfm.obnm,nx)
   if nx gt 0 then begin
      chiarrb[i]=cfb[i].mdchi
      chiarrm[i]=cfm[x].mdchi
      errarrb[i]=cfb[i].errvel
      errarrm[i]=cfm[x].errvel
   endif
endfor 

x=where(chiarrb ne 0,nx)
chiarrb=chiarrb[x]
chiarrm=chiarrm[x]
errarrb=errarrb[x]
errarrm=errarrm[x]

print,'number of observations:',nx
!p.charsize=1.8
!x.charsize=1.7
!y.charsize=1.7
plot,chiarrb,ps=8,col=1,xtitl='Obs number', ytit='chisq',/xsty
oplot,chiarrm,ps=8,col=222
!p.color=1
xyouts,100,1.0,/data,'!6 DSST from JJ conv code',size=1.9
!p.color=222
xyouts,100,0.8,/data,'!6 DSST from solar morph conv code',size=1.9
stop

plot,errarrb,ps=8,col=1,xtitl='Obs number', ytit='!6Errvel [ms!u-1!n]',/xsty
oplot,errarrm,ps=8,col=222
!p.color=1
xyouts,100,2.0,/data,'!6 DSST from JJ conv code',size=1.9
!p.color=222
xyouts,100,1.0,/data,'!6 DSST from solar morph conv code',size=1.9

end
