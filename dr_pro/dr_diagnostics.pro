pro dr_diagnostics, tag=tag, star=star,yrmo=yrmo

if ~keyword_set(star) then star = 'HR'
path='/tous/mir7/files_df/'

ff=file_search(path+'cd'+tag+'b'+star+'*'+yrmo+'*',count=count)
ncks=720
wavarr=fltarr(count,ncks)
delta_wav=fltarr(ncks, count) 

for i=0,count-1 do begin
   d1=strpos(ff[i],'cd')
   cbf=strmid(ff[i],d1,strlen(ff[i])-d1)
print,cbf
   dop_diagnostics,cdbfile=cbf,/wav_cont,/ctio,dwav=dwav
   delta_wav[*,i]=dwav
endfor 

print,'tag: ',tag
print,'median delta_wav: ',median(delta_wav)
print,'stddev delta_wav: ',stddev(delta_wav)
xlow=where(delta_wav lt 0.0,nxlow)
xhigh= where(delta_wav gt 2.0,nxhigh)
print,'number of cases where delta_wav lt 0.0: ',nxlow
print,'number of cases where delta_wav gt 2.0: ',nxhigh

stop
end
