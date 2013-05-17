pro tst_vank

for i=13,32 do begin
   vank,'128620','am','vdam',cf1,cf3,/ctio,ordr=[i,i]
endfor

num=32-13+1 
chi=fltarr(num)
err=fltarr(num)
vel=fltarr(num)

for i=0,num-1 do begin
   restore,'vstbank/vst128620_'+strcompress(string(i+13),/rem)+'.dat'
   chi[i]=cf3.mdchi
   err[i]=cf3.errvel
   vel[i]=cf3.mnvel
endfor

stop

end ;pro
