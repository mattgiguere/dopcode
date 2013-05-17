pro patch_d8, star=star

; counts correction for dewar 8
; output from vank is moved to cf2
; counts-corrected velocities go in cf1
; append prefix velocities to cf1 and store in cf3

if keyword_set(star) then ff=findfile('vstbank/vst'+star+'.dat',count=count) else $
    ff=findfile('vstbank/vst*.dat',count=count)

;counts fitting coefficients from decorel.pro
;mc=[9.47,  -0.0006,   3.19e-09]
;stop

for i=0,count-1 do begin
	restore,ff(i)
	cf2=cf1   ;save the unchanged dewar 8 mnvels in cf2
	; correct the cf1.mnvels
        for j=0,n_elements(cf1)-1 do begin
            if cf1(j).dewar eq 24 then begin
		tag=strmid(cf1(j).obnm,1,3)
		obnm=strmid(cf1(j).obnm,5,strlen(cf1(j).obnm)-5)
print,'tag: ',tag
print,'obnm: ',obnm
stop
;make a structure with 90th percentile counts for each observation

		path='/data/mir1_data/tag/'


	endif
        end
        cf3=cf1
        save,cf1,cf2,cf3,file=ff(i)
;	loadct,39
;        !p.background=255
;       velplot,cf1
;stop
end

;ck_pre,starnm=star

end
