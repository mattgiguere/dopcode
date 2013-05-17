pro dr_obs_lick_loop, star, obsnm=obsnm, pass=pass, vdavg=vdavg, debug=debug, demo=demo

; fischer 17nov2011
; dr_obs_lick_loop, '161797', obsnm='ri69.263', iss_nm='dsst161797fib_ri68.dat', tag='dff', 

  if keyword_set(demo) then demo=1 else demo=0

;	tag='dff'
;     	iss_nm='dsst161797fib_ri68.dat'   ;dsst created Jan 24, 2010
;     	iss_bc=-3831.739
;     	iss_obnm='ri68.109'
;       obsnm='ri69.'+strcompress(string(indgen(37)+86),/rem)   ;n1
;               obsnm='ri69.'+strcompress(string(indgen(50)+266),/rem)   ;n2a
;		vdavg='vddff_iod.ri69.n2a.avg'  ;n2a - implemented for pass3 only
;		vdavg='vddff_iod.ri69.n1.avg'   ;n1 - implemented for pass3 only
;               obsnm='ri68.'+strcompress(string(indgen(50)+48),/rem) 

;	tag='dfs'
;     	iss_nm='dsst161797slit_ri67.dat'   ;dsst created 15 Nov 2011
;     	iss_bc=-1100.229
;     	iss_obnm='ri67.287'
;            ;;;obsnm1='ri67.'+strcompress(string(indgen(48)+221),/rem)   
;            ;;;obsnm1='ri67.'+strcompress(string(indgen(42)+238),/rem)   
;       obsnm='ri70.'+strcompress(string(indgen(37)+76),/rem)   
;       vdavg='vddfs_iod.ri70.avg'   
;             obsnm=[obsnm1, obsnm2]

;	tag='dff'
;     	iss_nm='dsst188512fib_ri68.dat'   ;dsst created 15 Nov 2011
;     	iss_bc=11274.240
;     	iss_obnm='ri68.119'
;        obsnm='ri69.'+strcompress(string(indgen(33)+318),/rem)   
;        vdavg='vddff_iod.ri69.n2a.avg'   
;        vdavg='vddff_iod.ri69.n1.avg'   ;n1 - implemented for pass3 only

;	tag='dfs'
;     	iss_nm='dsst188512slit_ri67.dat'   ;dsst created Jan 24, 2010
;     	iss_bc=14702.777
;     	iss_obnm='ri67.297'
;        obsnm='ri70.'+strcompress(string(indgen(33)+128),/rem)   
;       vdavg='vddfs_iod.ri70.avg'   

;        tag='dfs'
;        iss_nm='dsst10700slit_rg92.dat'    ; dsst created 
;        iss_bc=4494.078
;        iss_obnm='rg92.225'
;        obsnm='ri70.'+strcompress(string(indgen(30)+218),/rem) 
;        obsnm='ri70.'+strcompress(string(indgen(20)+228),/rem) 
;        vdavg='vddfs_iod.ri70.avg'

        tag='dff'
        iss_nm='dsst10700fib_ri69.dat'    ; dsst created 
        iss_bc=18740.457
        iss_obnm='ri69.375'
;        obsnm='ri69.'+strcompress(string(indgen(32)+383),/rem) 
        obsnm='ri69.'+strcompress(string(indgen(32)+411),/rem) 
        vdavg='vddff_iod.ri69.n2c.avg'

num=n_elements(obsnm)

for i=0,num-1 do begin
     dopenv=lick_init(obsnm[i],iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag)
     restore,dopenv.bary_file   ;bcat
     x=where(bcat.obsnm eq obsnm[i],nx)
     if nx eq 0 then stop,'observation is not in qbcvel.dat' 
        found=findfile('/mir1/files/vd'+tag+star+'*'+obsnm[i], count=count)
        if count eq 1 then print, found
        if count eq 0 then begin
;           dop_main, star, dopenv, obsnm=obsnm[i], pass=1, tag=tag, $
;           		demo=demo;, verbose=verbose
;           dop_main, star, dopenv, obsnm=obsnm[i], pass=2, tag=tag
         dop_main, star, dopenv, obsnm=obsnm[i], pass=3, vdavg=vdavg, tag=tag
        endif
endfor

end
