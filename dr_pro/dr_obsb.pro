pro dr_obsb, run=run, obnm=obnm, pass=pass, debug=debug, demo=demo,tag=tag

star='128621'
if keyword_set(demo) then demo=1 else demo=0

  mpf=0  &   starsolve=0
  if ~keyword_set(lm) then lm='mpf'
  if lm eq 'mpf' then mpf=1
  if lm eq 'starsolve' then starsolve=1
  if ~keyword_set(lm) then mpf=1

  restore,'/mir7/bary/qbcvel.dat'
  xmatch=where(bcat.obsnm eq obnm,num)
  print,num
  if num eq 0 then stop

  obs_id=bcat[xmatch[0]].obsnm

     iss_nm='dsst128621m_rqa10.dat'
     iss_bc=-20540
     iss_obnm='rqa10.6622'
     dopenv=ctio_init(obs_id,iss_nm, iss_bc, iss_obnm=iss_obnm,tag=tag)
     tag='b'
     restore,dopenv.bary_file   ;bcat
     x=where(bcat.obsnm eq obs_id,nx)

    tag='b'
if nx gt 0 then begin
   bcat=bcat[x]
   for i=0,nx-1 do begin
      obnm = bcat[i].obsnm
      found=findfile('/mir7/files/vdb*'+obnm, count=count)
      if count eq 0 then begin
;         if keyword_set(debug) then $
;            dop_main, '128621',dopenv, obsnm=obnm, pass=1, tag=tag, /demo, /verbose else $
         dop_main, '128621', dopenv, obsnm=obnm,tag=tag,pass=1,demo=demo
         dop_main, '128621', dopenv, obsnm=obnm, tag=tag,pass=2
      endif
   endfor
endif


end
