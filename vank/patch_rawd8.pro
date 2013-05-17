pro patch_rawd8, star=star, percentile=percentile

; counts correction for dewar 8
; output from vank is moved to cf2
; counts-corrected velocities go in cf1
; append prefix velocities to cf1 and store in cf3

if keyword_set(star) then ff=findfile('vstbank/vst'+star+'.dat',count=count) else $
    ff=findfile('vstbank/vst*.dat',count=count)

restore,'comb_cts.dat'   ;cat

;counts fitting coefficients from decorel.pro
form1='(a8,i9,f9.2,i5)'

;stop
for i=0,count-1 do begin
	restore,ff(i)
	cf2=cf1   ;save the unchanged dewar 8 mnvels in cf2
	; correct the cf1.mnvels
        for j=0,n_elements(cf1)-1 do begin
            if cf1(j).dewar eq 24 then begin
               xm=where(cf1(j).obnm eq cat.obnm,nxm)
               if nxm eq 0 then begin 
                   print,cf1(j).obnm
                   mc=[9.47,  -0.0006,   3.19e-09]     
                    vel_cor=mc(0) + mc(1)*cf1(j).cts + mc(2)*cf1(j).cts^2  
                    if vel_cor gt 7.0 then vel_cor=7.  
                    if vel_cor lt -13.0 then vel_cor=-13.0    
                    print,vel_cor 
                endif
               if nxm gt 0 then begin
               if percentile eq 70 then begin
                   mc=[63.42, -0.03824]
                   bound=[1.4, -15.44]
                 vel_cor=mc(0) + mc(1)*cat(xm).ct_70 
 	         if vel_cor gt bound(0) then vel_cor=bound(0)
	         if vel_cor lt bound(1) then vel_cor=bound(1)
               endif
               if percentile eq 80 then begin
                   mc=[30.86, -0.0174]
                   bound=[2.44, -1752]
                 vel_cor=mc(0) + mc(1)*cat(xm).ct_80
 	         if vel_cor gt bound(0) then vel_cor=bound(0)
	         if vel_cor lt bound(1) then vel_cor=bound(1)
               endif
               if percentile eq 90 then begin
                   mc=[37.8449, -0.0227248, 2.3954e-06]
                   bound=[5.55, -15.963]
                 vel_cor=mc(0) + mc(1)*cat(xm).ct_90 + mc(2)*cat(xm).ct_90^2
 	         if vel_cor gt bound(0) then vel_cor=bound(0)
	         if vel_cor lt bound(1) then vel_cor=bound(1)
             endif
         endif
            print,cf1(j).obnm, cf1(j).cts, vel_cor, cf1(j).dewar,format=form1
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
