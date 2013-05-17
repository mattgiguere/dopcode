pro patch_d8, star=star

; counts correction for dewar 8
; output from vank is moved to cf2
; counts-corrected velocities go in cf1
; append prefix velocities to cf1 and store in cf3

if keyword_set(star) then ff=findfile('vstbank/vst'+star+'.dat',count=count) else $
    ff=findfile('vstbank/vst*.dat',count=count)

;counts fitting coefficients from decorel.pro
mc=[9.47,  -0.0006,   3.19e-09]
;stop
for i=0,count-1 do begin
	restore,ff(i)
	cf2=cf1   ;save the unchanged dewar 8 mnvels in cf2
	; correct the cf1.mnvels
        for j=0,n_elements(cf1)-1 do begin
            if cf1(j).dewar eq 24 then begin
               vel_cor=mc(0) + mc(1)*cf1(j).cts + mc(2)*cf1(j).cts^2
	       if vel_cor gt 7.0 then vel_cor=7.
	       if vel_cor lt -13.0 then vel_cor=-13.0
            print,vel_cor
            cf1(j).mnvel=cf1(j).mnvel-vel_cor
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
