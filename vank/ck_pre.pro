pro ck_pre,starnm=starnm

; tack on prefix velocities (SA accounted for by Chris)
; fischer aug2007

ff=findfile('vstbank/vst'+starnm+'.dat',count=count)
starname=strarr(count)

for i=0,count-1 do begin
	starname(i)=strmid(ff(i),11,strlen(ff(i))-15)
	restore,ff(i)
	n1=n_elements(cf1)
	n3=n_elements(cf3)
	if n1 eq n3 then begin
	  prenm='/mir1/vel/prebc/prebc'+starname(i)+'.dat'
	  preold='/mir1/vel/old.pre.includesSA/pre'+starname(i)+'.dat'
          pref=findfile(prenm,count=npre)
	  prefold=findfile(preold,count=npreold)
	  if npre eq 0 and npreold eq 1 then print,'no prebc file for: ',starname(i)
          if npre eq 1 then begin
	    print,'Prenosa file for: ',starname(i)
	    restore,pref
            xx=where(cf3(0).obnm eq pre.obnm,nxx)
	    if nxx eq 0 then if starname(i) eq 'gl285' $
	      or starname(i) eq '131156b' $
	      or starname(i) eq '95128' $
	      or starname(i) eq '9826' then begin 
		xx=where(cf3(1).obnm eq pre.obnm,nxx)
	 	diff =cf3(1).mnvel-pre(xx).mnvel
	    endif 
	    if nxx eq 1 then diff=cf3(0).mnvel-pre(xx).mnvel
            pre.mnvel=pre.mnvel+diff
            cf3=[pre,cf3(1:n_elements(cf1)-1)]
            ck1=n_elements(pre)+n_elements(cf1)-1
            ck2=n_elements(cf3)
            if ck1 ne ck2 then stop
	    save,cf1,cf2,cf3,f=ff(i)
	endif
	endif

endfor

end
