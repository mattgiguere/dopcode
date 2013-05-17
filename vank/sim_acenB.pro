pro sim_acenB

restore,'vst128621_cflat.dat'  ; linear trend from binary motion removed

prot_2011=36.71 ;d

e2s=3.0024584d-6  ; convert earth to solar masses
au=149597871000      ; AU in meters

t=cfc.jd        ; HJD
err=cfc.errvel
p=3.2357        ;d
tm=2455280.17   ;d
e=0.0
om=0.
k=0.51          ; m/s
a=0.04*AU       ; meters
inc=90.         ;deg
m1=0.934        ; solar mass units
m2=1.13*e2s     ; solar mass units
bigom=0.
cmvel=0.

rv, t, p, tm, e, om, a, inc, m1, m2, bigom, cmvel, vel1

;rv, t, prot_2011, tm, e, om, a, inc, m1, m2, bigom, cmvel, vel2
;rv, t, 6.2167, tm, e, om, a, inc, m1, m2*1.3, bigom, cmvel, vel3

jitter=sqrt(1.^2 + err^2)
 
;vel=vel1+vel2+vel3 
vel=vel1
dvel = randomn(seed, n_elements(vel))*jitter
vel=vel + dvel 
stop

plot, t, vel, ps=8, yra=[-10,10]
oploterr, t, vel, sqrt(0.6^2 + 0.7^2 + err^2)
stop

tim=t
data=vel
pergram,tim,data
plots,[p,p],[0,10],col=222,linesty=2,thick=4

stop

end ;pro 