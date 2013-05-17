pro mk_cf,starnm,tag=tag, obnm=obnm

if ~keyword_set(tag) then tag='f'
if ~keyword_set(starnm) then begin
 starnm=''  & tag=''
 read,'Enter the starname: ',starnm
 read,'Enter the tag for the cf file: ',tag
endif

restore,'reduce.dat'
x=where(cat.strnm eq starnm,nx)


if nx eq 0 and starnm ne 'MOON' then begin 
   addstar,star=starnm
   restore,'./reduce.dat'
   x=where(cat.strnm eq starnm,nx)
endif

   cf,starnm,cat(x).bccor,cfout
   cf=cfout
   nxx=0
   if keyword_set(obnm) then xx=where(strmid(cf.obnm,0,4) eq obnm,nxx)
   if nxx gt 0 then cf=cf[xx]
   save,cf,f='/mir7/dop/cf'+starnm+'_'+tag+'.dat'


end
