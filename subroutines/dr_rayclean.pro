pro dr_rayclean, starnm=starnm, date=date, savfile=savfile

;fischer aug2009

if ~keyword_set(starnm) then begin
   stnm=''
   read,'enter starname (e.g. 128620) ',stnm
   starnm=stnm
endif

if starnm eq '10700' then obsnm='achi120719.'+['1117','1118','1119']
if starnm eq '94683' then obsnm='achi120719.'+['1117','1118','1119']
if starnm eq '209100' then obsnm='achi120718.'+['1149','1150','1151']
if starnm eq '136442' then obsnm='achi120811.'+['1130','1131','1132']
if starnm eq '149661' then obsnm='achi120812.'+['1127','1128','1129']
if starnm eq '156274' then obsnm='achi120819.'+['1140','1141','1142']
if starnm eq '165341A' then obsnm='achi120820.'+['1136','1137','1138']
if starnm eq '22049' then obsnm='achi121103.'+['1141','1142','1143','1144','1145']
if starnm eq '20794' then obsnm='achi121104.'+['1127','1128','1129','1130','1131']


;128620
if starnm eq '128620' then begin
obsnm=['rqa10.6610','rqa10.6611','rqa10.6612','rqa10.6613','rqa10.6614',$
       'rqa10.6615','rqa10.6616','rqa10.6617','rqa10.6618','rqa10.6619',$
       'rqa10.6620','rqa10.6621','rqa10.6622','rqa10.6623','rqa10.6624',$
       'rqa10.6625','rqa10.6626','rqa10.6627','rqa10.6628','rqa10.6629',$
       'rqa10.6630','rqa10.6631','rqa10.6632','rqa10.6633','rqa10.6634']
endif

;128621
if starnm eq '128621' then begin 
obsnm=['rqa10.6645','rqa10.6646','rqa10.6647','rqa10.6648','rqa10.6649',$
       'rqa10.6650','rqa10.6651','rqa10.6652','rqa10.6653','rqa10.6654',$
       'rqa10.6655','rqa10.6656','rqa10.6657','rqa10.6658','rqa10.6659',$
       'rqa10.6660','rqa10.6661','rqa10.6662','rqa10.6663','rqa10.6664',$
       'rqa10.6665','rqa10.6666','rqa10.6667','rqa10.6668','rqa10.6669']
;       'rqa10.6670','rqa10.6671','rqa10.6672','rqa10.6673','rqa10.6674']
endif

if starnm eq '161797' then begin
obsnm=['ri67.280','ri67.281','ri67.282','ri67.283','ri67.284']

;obsnm=['ri68.104','ri68.105','ri68.106','ri68.107','ri68.108','ri68.109',$
;       'ri68.110','ri68.111','ri68.112','ri68.113']
endif

if starnm eq '188512' then begin
obsnm=['ri68.114','ri68.115','ri68.116','ri68.117','ri68.118','ri68.119',$
       'ri68.125','ri68.126','ri68.127','ri68.128']
endif

if ~keyword_set(savfile) then savfile='hd'+starnm+'.dat'

rayclean, obsnm, obstack, star, date=date,/auto


save,star,f='/tous/mir7/files_df/'+savfile

end


