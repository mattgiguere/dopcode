pro register

; nightly corrections based onshifts in common to A and B

restore,'/home/fischer/dop2/vank/vstbank/vst128620.dat'
cfa=cf3
x=where(cfa.jd gt 14800)  &  cfa=cfa[x]  ; no B obs run before this
;remove linear trend
coefa=poly_fit(cfa.jd,cfa.mnvel,1)  & cfa.mnvel=cfa.mnvel-poly(cfa.jd,coefa)
velplot,cfa,titl='aCenA',9.5/24.,d,e,v,n,c,duma  & cora=duma
numa=n_elements(duma)
jd=fltarr(numa)  &  cor_ab=fltarr(numa)

restore,'/home/fischer/dop2/vank/vstbank/vst128621.dat'
cfb=cf3
x=where(cfb.jd gt 14800)  &  cfb=cfb[x]  ; no B obs run before this
;remove linear trend
coefb=poly_fit(cfb.jd,cfb.mnvel,1)  & cfb.mnvel=cfb.mnvel-poly(cfb.jd,coefb)
velplot,cfb,titl='acenB',9.5/24.,d,e,v,n,c,dumb  & corb=dumb

;stop
for i=0,numa-1 do begin
   jd[i]=round(duma[i].jd)
   xa=where(round(duma.jd) eq jd[i],nxa)
;   if nxa ne 1 then stop
   xb=where(round(dumb.jd) eq jd[i],nxb) 
   if nxb gt 0 then begin   ;average offset for AB
      cor_ab[i]=mean([duma[xa].mnvel,dumb[xb].mnvel])
      cora[xa].mnvel=cora[xa].mnvel-cor_ab[i]
      corb[xb].mnvel=corb[xb].mnvel-cor_ab[i]
   endif
;   if nxa gt 0 and nxb eq 0 then begin  ; only offset for A
;      cor_ab[i]=duma[xa].mnvel
;      cora[xa].mnvel=cora[xa].mnvel-cor_ab[i]
;   endif
;   if nxb gt 0 and nxa eq 0 then begin  ; only offset for A
;      cor_ab[i]=dumb[xb].mnvel
;      cora[xa].mnvel=cora[xa].mnvel-cor_ab[i]
;   endif
endfor 
;stop
plot,duma.jd,duma.mnvel,ps=8,col=0
oplot,dumb.jd,dumb.mnvel,ps=8,col=222
;stop

plot,cora.jd,cora.mnvel,ps=8,col=0
oplot,corb.jd,corb.mnvel,ps=8,col=222

;correct vst for A
restore,'vstbank/vst128620.dat'  &  cfa=cf3
xaa=where(cfa.jd gt 14800.)  & cfa=cfa[xaa]
for i=0,n_elements(cfa.jd)-1 do begin
   xxa=where(round(cfa[i].jd) eq jd,nxxa)
   if nxxa gt 0 then cfa[i].mnvel = cfa[i].mnvel-cor_ab[xxa[0]]
endfor

;correct vst for B
restore,'vstbank/vst128621.dat'  &  cfb=cf3
xbb=where(cfb.jd gt 14800.)  & cfb=cfb[xbb]
for i=0,n_elements(cfb.jd)-1 do begin
   xxb=where(round(cfb[i].jd) eq jd,nxxb)
   if nxxb gt 0 then cfb[i].mnvel = cfb[i].mnvel-cor_ab[xxb[0]]
endfor

velplot,cfa


stop
end
