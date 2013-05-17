pro ck_wdif, star=star, tag=tag, yrmo=yrmo, all=all

if keyword_set(yrmo) then ff=file_search('/tous/mir7/files_df/cd'+tag+'b'+star+'_achi'+yrmo+'*',count=count)
if keyword_set(all) then ff=file_search('/tous/mir7/files_df/cd'+tag+'b'+star+'_achi*',count=count)

med_wdif0=fltarr(count)
med_ddif0=fltarr(count)
med_wdif1=fltarr(count)
med_ddif1=fltarr(count)
obnm=strarr(count)
num=760  ; n_chunks

	wdif0=fltarr(num)
	ddif0=fltarr(num)
	wdif1=fltarr(num)
	ddif1=fltarr(num)

for i=0,count-1 do begin
	restore,ff[i]
	d1=strpos(ff[i],'achi')
	obnm[i]=strmid(ff[i],d1,strlen(ff[i])-d1)
	for j=0,num-1 do begin
		wdif0[j]=(*chunk_arr[j].free_par)[0].wcof[0] - (*chunk_arr[j].free_par)[2].wcof[0]  
		ddif0[j]=(*chunk_arr[j].free_par)[0].wcof[1] - (*chunk_arr[j].free_par)[2].wcof[1]  
		wdif1[j]=(*chunk_arr[j].free_par)[1].wcof[0] - (*chunk_arr[j].free_par)[2].wcof[0]  
		ddif1[j]=(*chunk_arr[j].free_par)[1].wcof[1] - (*chunk_arr[j].free_par)[2].wcof[1]  
	endfor
	med_wdif0[i]=median(wdif0)
	med_ddif0[i]=median(ddif0)
	med_wdif1[i]=median(wdif1)
	med_ddif1[i]=median(ddif1)
endfor 

!x.margin=[2,2]
!y.margin=[2,2]
if strmid(star,0,2) eq 'HR' then begin
	!p.multi=[0,2,2]
	plot,med_wdif0,xtitl='wav diff 0',ps=8
	plot,med_wdif1,xtitl='wav diff 1',ps=8
	plot,med_ddif0,xtitl='disp diff 0',ps=8
	plot,med_ddif1,xtitl='disp diff 1',ps=8
endif
!p.multi=[0,1,1]

if strmid(star,0,2) ne 'HR' then begin
	flag=intarr(count)
	vel=fltarr(count)
	restore,'/Users/debrafischer/dop2/vank/vstbank/vst'+star+'.dat'
	for i=0,count-1 do begin
		x=where(cf3.obnm eq obnm[i],nx)
		if nx eq 1 then flag[i]=1
		if nx eq 1 then vel[i]=cf3[x].mnvel
	endfor
	!p.multi=[0,2,2]
	xgd=where(flag eq 1)
	plot,med_wdif0[xgd],xtitl='wav diff 0', vel[xgd],ps=8
	plot,med_wdif1[xgd],xtitl='wav diff 1',vel[xgd],ps=8
	plot,med_ddif0[xgd],xtitl='disp diff 0', vel[xgd],ps=8
	plot,med_ddif1[xgd],xtitl='disp diff 1', vel[xgd],ps=8
endif
!p.multi=[0,1,1]

print,'median change in wave for '+star, median(med_wdif0),median(med_wdif1)
print,'median change in disp for '+star, median(med_ddif0),median(med_ddif1)
	
stop
end  ;pro 

