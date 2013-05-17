pro rv,t,P,tp,e,om,a,i,m1,m2,bigom,cmvel,vel
;
; NAME:
;
;  RV
;
; PURPOSE:
;
;  This procedure will return the heliocentric radial velocities of a
;  spectroscopic double lined binary system at a given heliocentric Julian
;  date and specified orbital parameters.
;
; CALLING SEQUENCE:
;
;  rv,t,P,tp,e,om,a,i,m1,m2,bigom,cmvel,vel
;
; INPUT:
;
;  t	-time at each position of orbit calculated (scalar or vector)
;	[reduced HJD]
;  P	-period [days]
;  tp	-time of last periapse [reduced HJD]
;  e 	-eccentricity of the orbit
;  om	-line of apsides [deg]
;  a	-semi-major axis (relative) [m]
;  i	-angle between the line of sight and the axis of orbit [deg]
;  m1	-mass of primary star [solar masses]
;  m2	-mass of secondary star [solar masses]
;  bigom-longitude of ascending node [deg]
;  cmvel-center of mass velocity [m/s]
;
; OUTPUT:
;
;  VEL  fltarr  Predicted heliocentric radial velocities [m/s] 
;               for the date(s) specified.
;  max radial velocities 	-max rv1 and rv2, m/s
;  closest approach 		-smallest seppararion, solar radii
;  vel at closest approach	-rv2 at closest approach, m/s
;  time of closest approach	-reduced HJD of closest approach
;
;  
; RESTRICTIONS:
;  
;  To ensure consistency with the routines JULDATE and HELIO_JD, the 
;  reduced HJD must be used throughtout.
;
; EXAMPLES:
; 
;
;
; 
;
; Check that a suitable number of parameters have been entered.
   if n_params() lt 1 then begin
    print,'syntax: rv,t,P,tp,e,om,a,i,m1,m2,bigom,cmvel,vel
    retall
   endif

; Check to make sure that reduced HJD has been used 
;   if (max(t) ge 1.d6) then message,'Full HJD entered, use reduced HJD.'
;   if (tp ge 1.d6) then message,'Full HJD entered, use reduced HJD'

; Number of positions to be calculated
   nt=n_elements(t)

; Constants:
  AU=1.49597870662d11			;m
  G=6.6725985d-11			;m^3/(kg*s^2)
  Msun=1.9889225d30			;kg
  day=24.d0 * 60.d0 * 60.d0		;seconds, dbl precision

; Convert arguments to double precision.

   P=double(P)		;period, days
   tp=double(tp)
   a=double(a)		;semimajor axis (relative), m
   e=double(e)	;eccentricity
   i=double(i*!dtor)	;inclination, radians
   om=double(om*!dtor)		;longitude of ascending nodes, radians
   bigom=double(bigom*!dtor)	;line od apsides, radians
   cmvel=double(cmvel)		;center of mass velocity, m/s
   m1=double(max([m1,m2]))	;mass1, solar masses
   m2=double(min([m1,m2]))	;mass2, solar masses
   mt=double(m1+m2)		;total mass, solar masses
   if n_elements(palallax) ne 1 then palallax=0

; Find the missing parameter.
   if mt le 0 then begin			;no mass function
    mt=((2.d0*!dpi)/P)^2.d0*a^3/G	
    print,'Sum of masses (solar masses):	',mt/Msun	
   endif
   if P le 0 then begin				;no period
    P=2.d0*!dpi*a^1.5/(G*mt)^0.5
    print,'Orbital period (days):		',P/day
   endif
   if a le 0 then begin				;no semimajor axis (rel)
    a=(G*mt*(P/2.d0*!dpi)^2.d0)^(1.d0/3.d0)
    print,'Semimajor axis (rel) (m)		',a/m
   endif

; Calculate the approximate eccentric anomaly, E1, via the mean
; anomaly, M.
   M=2.d0*!dpi*( ((t-tp)/P) - fix((t-tp)/P))
;   E1=M+e*sin(M)+((e^2)*sin(2.d0*M)/2.d0)+((e^3)*3.d0^2*sin(3.d0*M) $
;      -3.d0*sin(M)/(24.d0))
;eorig=e1
;New First Guess: Taff pg.54
   E1 = M + e*sin(M)/(1.d0-sin(M+e)+sin(M))

;  Ea=E0-e*sin(E0)-M
;  Eb=1.d0-e*cos(E0)

; Refine the estimate using the Newton-Rhapson iteration
   ct = 0
   repeat begin
    ct = ct+1
    E0=E1
    Ea=E0-e*sin(E0)-M
    Eb=1.d0-e*cos(E0)
    E1=E0+(Ea/Eb)*(E1 gt E0)-(Ea/Eb)*(E1 le E0) 
;    if ct gt 100 then begin
;      dum = max(E1-E0,i)
;      print,'E1=',E1
;    end
   endrep until (max(E1-E0) lt 1.0d-6 or ct gt 100)

; Calculate nu and r
   n1=1.d0 + e
   n2=1.d0 - e
   nu=2.d0*atan((n1/n2)^0.5*tan(E1/2.d0))
   r=a*(1.d0-e*cos(E1))				;in meters

; Calculate the radial velocity
   k1=(2.d0*!dpi/P)*(m2*a/(m1+m2))*sin(i)/(1.d0-(e^2))^0.5
   rv1=k1*(cos(nu+om)+e*cos(om))
;   k2=(2.d0*!dpi/P)*(m1*a/(m1+m2))*sin(i)/(1.d0-(e^2))^0.5
;  rv2=k2*(cos(nu+om)+e*cos(om))

; convert the radial velocity from meters/day to meters/sec
   vel=rv1/day+cmvel
;   rv2=rv2*AU/day

; plot the result
;   if keyword_set (plot) then begin
;    ymn=min([rv1,rv2],max=ymx)
;    xtit='Time (days)'
;    ytit='Radial Velocity (m/s)'
;    plot,to,rv1,xtit=xtit,ytit=ytit,xsty=2,yr=[ymn,ymx],/nodata
;    oplot,to,rv1,co=6
;    oplot,to,rv2,co=2
;    oplot,!x.crange,[0,0],line=3,co=3
;   endif

;  max1=max(rv1)
;  max2=max(rv2)
;  print,'Max radial velocities (m/s):		',max1>max2, max1<max2
;  print,'Closest approach (solar radii):	',min(r,min)/6.96d8
;  print,'Velocity at closest approach (m/s):	',rv2(min)
;  print,'Time of closest approach (days):	',to(min)/day

end
   







