pro ck_vel

restore,'vst75732_da111706_A.dat'
cfd=cf1
num=n_elements(cfd)

diff_rv=fltarr(num)
diff_err=fltarr(num)

;restore,'/mir1/vel/vst75732.dat'
restore,'vst75732.dat'

for i=0,num-1 do begin 
	x=where(cfd(i).obnm eq cf1.obnm,nx)
	if nx eq 0 then print,cfd(i).obnm
	diff_rv(i)=cfd(i).mnvel-cf1(x).mnvel
	diff_err(i)=cfd(i).errvel-cf1(x).errvel
end

plothist,diff_rv

stop
end








