pro dr_keck

readcol,'starnames_kecktempl.txt',name,f='a9'
num=n_elements(name)

restore,'/mir1/lick_st.dat'

for i=0,num-1 do begin
    vank,name(i),'f','vdf',mct=200,fit_key=2
print,name(i)
x=where(name(i) eq lick.name,nx) 
if nx eq 0 then print,name(i)+' not in lick structure'
stop
end



end
