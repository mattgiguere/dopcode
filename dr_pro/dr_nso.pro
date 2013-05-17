pro dr_nso, starnm, run=run, obnm=obnm, tag=tag, pass=pass, debug=debug

restore,'/data/ctio/bary/qbcvel.dat'
if keyword_set(run) then x=where(strmid(bcat.obsnm,0,5) eq run $
                                 and bcat.objnm eq starnm and $
                                 bcat.obtype eq 'o',nx)
if keyword_set(obnm) then x=where(bcat.obsnm eq obnm and $
                                   bcat.obtype eq 'o',nx)
tag='b'

if nx gt 0 then begin
   bcat=bcat[x]
   for i=0,nx-1 do begin
      obnm = bcat[i].obsnm
      found=findfile('/data/ctio/files/vd'+tag+'*'+obnm, count=count)
      if count eq 0 then begin
         if keyword_set(debug) then $
            dop_main, starnm, obsnm=obnm, pass=1, /demo,/nso else $
         dop_main, starnm, obsnm=obnm, pass=1, tag=tag, /nso
         dop_main, starnm, obsnm=obnm, pass=2, tag=tag, /nso

      endif
   endfor
endif


end
