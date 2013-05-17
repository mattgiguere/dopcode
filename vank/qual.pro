pro qual

newf=findfile('vst*.dat',count=count)
discrep=''

dvel=fltarr(count)
derr=fltarr(count)
       
;!p.multi=[0,3,1]
loadct,39
!p.background=255

for i=0,count-1 do begin
    oldf=findfile('/mir1/vel/'+newf(i),count=ct)
    if ct eq 1 then begin
        starnm=strmid(newf(i),3,strlen(newf(i))-7)
        print,starnm
        restore,oldf
        cfold=cf1
        restore,newf(i)
        cfnew=cf1
        num=n_elements(cfnew)
        dv=fltarr(num)  & de=fltarr(num) & ob=strarr(num) & chi=fltarr(num)
        for ii=0,num-1 do begin
            dv(ii)=9999.
            de(ii)=9999.
            tt=where(cfnew(ii).obnm eq cfold.obnm,ntt)
            if ntt eq 1 then begin 
               dv(ii)=cfnew(ii).mnvel-cfold(tt).mnvel
               de(ii)=cfnew(ii).errvel-cfold(tt).errvel
               ob(ii)=cfnew(ii).obnm
               chi(ii)=cfnew(ii).mdchi-cfold(tt).mdchi
           endif
       endfor
       xg=where(dv ne 9999.,nxg)
       if nxg gt 0 then begin
           dv=dv(xg)
           de=de(xg)
           ob=ob(xg)
           chi=chi(xg)
       endif
       save,dv,de,ob,filename=starnm+'_ck.dat'
;!p.color=1
;       plothist,dv,bin=.2,yra=[0,5]
;       plothist,de,bin=.2,yra=[0,5]
;       plothist,chi,bin=0.2,yra=[0,5]
       hlp=where(abs(dv) gt 3,nhlp)
       if nhlp gt 0 then begin
           print,'Discrepant points: ',ob(hlp)
           discrep=[discrep,ob(hlp)]
       endif

   end
end

stop
save,discrep,'discrep.dat'

end

