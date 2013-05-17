pro ckiodt,fts=fts,tag=tag

; checks the chisq fit for a series of test iodine observations 
; taken with different temperature settings on the temp controller


if ~keyword_set(tag) then read,'enter the run tag (b): ',tag
if fts eq 'eso50' then r1='vd'+tag+'iod_IODINE_'
if fts eq 'pnnl' then r1='vd'+tag+'iod_IODINE_'
r2='_rqa10.'

t35=r1+'35C'+r2+['7151','7152','7153','7154','7155','7156',$
              '7157','7158','7159','7160']

t40=r1+'40C'+r2+['6801','6802','6803','6804','6805','6806',$
             '6807','6808','6809','6810']

t45=r1+'45C'+r2+['6811','6812','6813','6814','6815','6816',$
             '6817','6818','6819','6820']

t50=r1+'50C'+r2+['6821','6822','6823','6824','6825','6826',$
             '6827','6828','6829','6830']

t55=r1+'55C'+r2+['6831','6832','6833','6834','6835','6836',$
             '6837','6838','6839','6840']

t60=r1+'60C'+r2+['6841','6842','6843','6844','6845','6846',$
             '6847','6848','6849','6850']

t65=r1+'65C'+r2+['6900','6901','6902','6903','6904','6905',$
             '6906','6907','6908','6909']

t70=r1+'70C'+r2+['6910','6911','6912','6913','6914','6915',$
             '6916','6917','6918','6919']

path='/data/ctio/files/'
chi_35=fltarr(10)   &  chi_40=fltarr(10)  &  chi_45=fltarr(10) 
chi_50=fltarr(10)   &  chi_55=fltarr(10)  &  chi_60=fltarr(10) 
chi_65=fltarr(10)   &  chi_70=fltarr(10)  

for i=0,9 do begin
   ff=findfile(path+t35[i],count=count)
   if count gt 0 then restore,ff
   chi_35[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t40[i],count=count)
   if count gt 0 then restore,ff
   chi_40[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t45[i],count=count)
   if count gt 0 then restore,ff
   chi_45[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t50[i],count=count)
   if count gt 0 then restore,ff
   chi_50[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t55[i],count=count)
   if count gt 0 then restore,ff
   chi_55[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t60[i],count=count)
   if count gt 0 then restore,ff
   chi_60[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t65[i],count=count)
   if count gt 0 then restore,ff
   chi_65[i]=median(vd.fit)
endfor

for i=0,9 do begin
   ff=findfile(path+t70[i],count=count)
   if count gt 0 then restore,ff
   chi_70[i]=median(vd.fit)
endfor

x1=where(chi_35 ne 0.0)  & chi_35=chi_35[x1]
x2=where(chi_40 ne 0.0)  & chi_40=chi_40[x2]
x3=where(chi_45 ne 0.0)  & chi_45=chi_45[x3]
x4=where(chi_50 ne 0.0)  & chi_50=chi_50[x4]
x5=where(chi_55 ne 0.0)  & chi_55=chi_55[x5]
x6=where(chi_60 ne 0.0)  & chi_60=chi_60[x6]
x7=where(chi_65 ne 0.0)  & chi_65=chi_65[x7]
x8=where(chi_70 ne 0.0)  & chi_70=chi_70[x8]

mc35=mean(chi_35)  &  print,'Mean chi T=35: ',mc35
mc40=mean(chi_40)  &  print,'Mean chi T=40: ',mc40
mc45=mean(chi_45)  &  print,'Mean chi T=45: ',mc45
mc50=mean(chi_50)  &  print,'Mean chi T=50: ',mc50
mc55=mean(chi_55)  &  print,'Mean chi T=55: ',mc55
mc60=mean(chi_60)  &  print,'Mean chi T=60: ',mc60
mc65=mean(chi_65)  &  print,'Mean chi T=65: ',mc65
mc70=mean(chi_70)  &  print,'Mean chi T=70: ',mc70


stop   
end
