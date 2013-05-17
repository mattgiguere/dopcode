pro dr_iod_test, demo=demo, tag=tag, run=run

if ~keyword_set(demo) then demo=0

;anss=''
;read,'Are you sure the psf description is properly set? (y/n) ',anss
;if anss eq 'n' then stop

;obsnm_arr='rj162.'+['63','64','65','81','82','83','84','85','86','92','93','94',$
;	'139','140','141','147','148','149','150','176','177','178','184','185','186',$
;	'217','218','219','225','226','227','243','244','245','251','252','253',$
;	'269','270','271','277','278','279','285','286','287']

;obsnm_arr='rj162.'+['350','351','352','403','404','405','430','431','432']

restore,'/tous/mir3/bary/kbcvel.dat'
x=where(strmid(bcat.obsnm,0,5) eq run and strmid(bcat.objnm, 0,2) eq 'HR') 

obsnm_arr=bcat[x].obsnm 
stop
num=n_elements(obsnm_arr) 

for i=0,num-1 do begin	
	dr_iod_soln,obsnm_arr[i],tag='js',demo=demo,/keck
endfor
 
end
