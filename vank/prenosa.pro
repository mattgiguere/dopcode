pro prenosa,star=star

restore,'/home/fischer/vank/vstbank/vst'+star+'.dat'
restore,'/mir1/vel/prebc/prebc'+star+'.dat'

num=n_elements(cf1)

x=where(cf3(0).obnm eq pre.obnm,nx)
if nx eq 0 then stop

if nx eq 1 then begin
	diff=cf3(0).mnvel-pre(x).mnvel
	pre.mnvel=pre.mnvel+diff
	cf3=[pre,cf3(1:num-1)]
	ck1=n_elements(pre)+n_elements(cf1)-1
	ck2=n_elements(cf3)
	if ck1 ne ck2 then stop
end

save,cf1,cf3,f='/home/fischer/vank/vstbank/vst'+star+'.dat'

end
