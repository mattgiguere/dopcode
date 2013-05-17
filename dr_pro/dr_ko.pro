pro dr_ko, star=star

restore,'/tous/mir3/bary/kbcvel.dat'
x=where(strmid(bcat.obsnm,0,5) eq 'rj162' and bcat.objnm eq star and $
  bcat.obtype eq 'o',nx)
print,'the following I2 observations were found: '
print,bcat[x].obsnm 
obs_obnm=bcat[x].obsnm

for i=0,nx-1 do begin
   dr_obskeck,obnm=obs_obnm[i],star=star
endfor


stop

end  ;pro
