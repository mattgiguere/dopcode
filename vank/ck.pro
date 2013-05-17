pro ck

restore,'vst75732_daf.dat'
cfd=cf1
restore,'vst75732.dat
cfo=cf1

numd=n_elements(cfd)
numo=n_elements(cfo)

for i=0,numd-1 do begin
	x=where(cfd.obnm eq cfo.obnm,nx)
	if nx eq 0 then print,cfd(i).obnm
end

stop
end
