pro ck_hip8102

restore,'vstbank/vsthip8102.dat'

num=n_elements(cf1)

for i=0,19 do begin
    plot,cf1.mdpar(i),cf1.mnvel,ps=8
    wait,2
end

    plot,cf1.mdchi,cf1.mnvel,ps=8,xtit='Chi', ytit='mnvel'



;for i=0,num-1 do begin
;    plot,cf1.mnpar(i),cf1.mnvel,ps=8
;    stop
;end

stop

end

