
; Bread and butter weight routine.
function assign_weights5, vdarr, nchunk, numobs, fitdex, diff, sigma, $
	max_chi=max_chi, cut_thres=cut_thres
  ;fischer nov09
  ; 1. want weight to decrease as chisq fit increases 
  ;      - toss chunk sets with bad chisq fits (telluric lines) 
  ; 2. want weight to decrease as difference from the median vel in a chunkset increases
  ;Find scatter of diff within a chunk set
  ;The ``sigma'' below measures the SCATTER of chunk set members relative
  ;to their expected velocity based on the observation median.

  fudgefac = 1.      ;double the sigmas for each chunk
;  fudgefac = 2.      ;double the sigmas for each chunk
  tmp_fit=(vdarr.(fitdex))
  i=sort(tmp_fit)
  tmp_fit=tmp_fit[i]
  if cut_thres lt 0.6 or cut_thres gt 0.9 then cut_thres=0.95  ;need sensible value
  findx=[cut_thres*n_elements(tmp_fit)]
  fitlimit=tmp_fit[findx]
  fitlimit=fitlimit[0]   ; real number
  print,'fit limit: ',fitlimit
  sig_check=stddev(diff)  
  print,'unaltered stddev: ', sig_check
  if ~keyword_set(max_chi) then max_chi=3.
  FOR n = 0,nchunk-1 do begin
     igood = where(vdarr(*,n).(fitdex) lt max_chi) ;obs having good fit
     sigma(n) = stddev(diff(igood,n))                     ;RMS of vel-median_obsof chunk-set 
  END

  FOR ob = 0,numobs-1 do begin
      const = median( abs(diff(ob,*))/sigma)   ;ratio: actual to typical discrepancy for chunks
      sigmaob = const * sigma * fudgefac       ;boost chunk-set sigma by const
      for n=0,nchunk-1 do begin                ;hi chisq fits drop out
        vdarr(ob,n).weight = 1./sigmaob(n)^2   ;weights scaled to scatter w/i obs.
      end
      ;vdarr(*,n).weight = 1.             ;No weights
      ; Otherwise, weights based on input (hopefully photon limited errors)
  END

  return, 0
end; assign_weights5

function assign_weights1, vdarr, nchunk, numobs, fitdex, diff, sigma, max_chi=max_chi
  ;Find scatter of diff within a chunk set
  ;The ``sigma'' below measures the SCATTER of chunk set members relative
  ;to their expected velocity based on the observation median.
  if ~keyword_set(max_chi) then max_chi=3.
  fudgefac = 1.      ;double the sigmas for each chunk
;  fudgefac = 2.      ;double the sigmas for each chunk
  FOR n = 0,nchunk-1 do begin
      igood = where(vdarr(*,n).(fitdex) lt max_chi)  ;obs having fit better than 3.
      sigma(n) = stddev(diff(igood,n))         ;RMS of vel-median_obsof chunk-set 
;      ibad=where(vdarr(*,n).(fitdex) gt 10.,nibad)
;      if nibad gt 0 then stop
  END
 
  FOR ob = 0,numobs-1 do begin
      const = median( abs(diff(ob,*))/sigma )  ;ratio: actual to typical discrepancy for chunks
      sigmaob = const * sigma * fudgefac       ;boost chunk-set sigma by const
      for n=0,nchunk-1 do begin
        vdarr(ob,n).weight = 1./sigmaob(n)^2   ;weights scaled to scatter w/i obs.
        if vdarr[ob,n].fit gt max_chi then vdarr[ob,n].weight = 0.0d
      end
      ;vdarr(*,n).weight = 1.             ;No weights
      ; Otherwise, weights based on input (hopefully photon limited errors)
  END

  return, 0
end; assign_weights1

; Weighting routine tried on 2/1/06 - dave
function assign_weights2, vdarr, diff,  nchunk, fitdex, numobs
  ; dave test - uncomment above when done
  print, 'DAVE TEST DAVE TEST DAVE TEST DAVE TEST'
  for n=0,nchunk-1 do begin
      for ob=0, numobs-1 do begin
          igood1 = where(vdarr[*,n].(fitdex) lt 99)

          ; remove this observation
          q = where(igood1 eq ob)
          if (q ne -1) then begin
              first = 0
              last = n_elements(igood1)-1
              CASE q OF
                  first: igood2 = igood1[1:n_elements(igood1)-1]
                  last: igood2 = igood1[first:last-1]
                  ELSE: igood2 = [ igood1[first:q-1], igood1[q+1:last] ]
              ENDCASE
          endif

          ; Calculate the stddev of the other observations in this set.
          sig = stddev(diff[igood2,n])

          vdarr[ob,n].weight = (1./sig^2)

          ; Reduce the weight if the diff for this chunk is too far away from the majority.
          if (abs(diff[ob,n]) gt sig) then begin
              factor = (abs(diff[ob,n]) / sig)^2
              vdarr[ob,n].weight = vdarr[ob,n].weight * (1./factor)
          endif
      endfor

      ;window, 1
      ;plothist, diff[*,n], bin=10, xr=[-150,150]
      ;str = 'chunk set ' + strcompress(string(n))
      ;xyouts, 10, 10, str

      ;window, 2
      ;plot, vdarr[*,n].weight, psym=4

      ;print, 'stdev = ', stddev(diff[igood1,n])

      ;stop
  endfor

  return, 0
end

; Weighting scheme from 2/8/06 - dave
function assign_weights3, vdarr, nchunk, fitdex, numobs 
  for ob=0,numobs-1 do begin
      minfit = min(vdarr[ob,*].(fitdex))
      maxcount = max(vdarr[ob,*].cts)
      for n=0,nchunk-1 do begin
          chunk = vdarr[ob,n]
          vdarr[ob,n].weight = (minfit / chunk.(fitdex))^2 * (double(chunk.cts) / maxcount)
      ;    vdarr[ob,n].weight = (minfit / chunk.(fitdex))
      endfor

  endfor
end

; Weighting scheme from 2/10/06 - dave
function assign_weights4, vdarr, nchunk, numobs, fitdex
  for ob=0,numobs-1 do begin
      for n=0,nchunk-1 do begin
          chunk = vdarr[ob,n]
          vdarr[ob,n].weight = 1 / abs(vdarr[ob,n].ivel - 0)
      endfor

  endfor
end


pro vdclean, vdarr,ifit,nchunk,plot=plot,maxchi=maxchi, nozero=nozero, sp2=sp2
;
;  Cleans an array of VD structures, VDCUBE, of "bad" chunk sets.
;  "Bad" chunks are those having:
;   1.  high photon-limited errors
;   2.  high ChiSq, on average over all chunks
;  The underlying idea is that some chunks are bad, due to
;  problems with the DSST (==> chi poor) or poor Doppler information 
;  (small  slopeweight).
;
;INPUT:
;  vdarr    (output structure)  contains an array of VD structures, for 
;                               each observation and chunk within it.
;  fitdex   (integer)           the index of the chisq fit (9 or 10 or ...)
;
;OUTPUT:
; vdarr     (output structure)  contains a "cleaned" array of VD structures .
; nchunk                        revised (reduced) # of chunks
;
;OPTIONAL:
; plot      (keyword)           plots histogram of "goodness" parameter
;

;res=finite(vdarr[0,*].vel) 
;gg=where(res eq 0,ngg)
;if ngg gt 0 then print, 'NaN: ',gg

IF n_params () lt 2 then begin
  print,'-------------------------------------------------------------------'
  print,' SYNTAX:'
  print,' '
  print,' IDL> vdclean, vdarr,ifit,nchunk,plot=plot'
  print,'-----------------------------------------------------------------'
  RETURN
ENDIF


;   
;Toss sets with vd.weight lt 0, Bad DSST chunks, PB 6/17/98
print,'VDCLEAN, input n_chunks '+string(n_elements(vdarr(0,*)))
igood = where(vdarr[0,*].weight gt 0)
vdarr = vdarr(*,igood)      ;use good chunk sets only
print,'VDCLEAN, output n_chunks '+string(n_elements(vdarr(0,*)))

;SET CONSTANTS and DEFAULTS
     tagnam = tag_names(vdarr)
     orignum = n_elements(vdarr(0,*).vel)  ;# of chunks in observation #0
     numobs = n_elements(vdarr(*,0))        ;# of observations
     mdchi = fltarr(orignum)
     if not keyword_set(maxchi) then maxchi = 10.
     fitdex = first_el(where(tagnam eq 'FIT')) 
     veldex = first_el(where(tagnam eq 'VEL')) 
     if ifit eq 1 then begin
        fitdex = first_el(where(tagnam eq 'FIT'))
        veldex = first_el(where(tagnam eq 'VEL'))
        if veldex lt 0 then veldex = first_el(where(tagnam eq 'VEL'))
     endif
     if ifit eq 2 then begin
        fitdex = first_el(where(tagnam eq 'SFIT'))
        veldex = first_el(where(tagnam eq 'SVEL'))
        if veldex lt 0 then veldex = first_el(where(tagnam eq 'VEL'))
     endif
     if ifit eq 3 then begin
        fitdex = first_el(where(tagnam eq 'FIT1'))
        veldex = first_el(where(tagnam eq 'VEL1'))
        if veldex lt 0 then veldex = first_el(where(tagnam eq 'VEL'))
     endif

;Reject chunks having: high photon-limited error, high CHISQ, or both

err = 1./sqrt(vdarr(0,*).weight)         ;photon-limited error (m/s)

FOR ch = 0, orignum-1 do begin
   mdchi(ch) = median(vdarr(*,ch).(fitdex)) ;median(chisq) 
   if mdchi(ch) eq 0 then begin
     print,'Observation has fit/ifit/sfit = 0'
     print,'Observation probably has not sucessfully made a second pass' 
;    mdchi(ch)  = median(vdarr(*,ch).fit)  ;insert when one-pass desired
     return
   endif
endfor
err=reform(err)
ierr = sort(err)   &  therr  = err(ierr(0.96*orignum))      ;99 %'ile.
ichi = sort(mdchi) &  thchi = mdchi(ichi(0.96*orignum))      ;99 %'ile
;
if keyword_set(maxchi) then thchi = min([maxchi,thchi])

;NEW WEIGHTS, include median of CHISQ for each chunk-set.
;for ob=0,numobs-1 do vdarr(ob,*).weight = 1./(err^2*mdchi)  ;new weights
;print,'WEIGHTS=1' & for ob=0,numobs-1 do vdarr(ob,*).weight = 1.  ;new weights

igood = where(err le therr   and $          ;accept low photon errors
              mdchi lt thchi, nchunk)               ;accept low CHISQ

vdarr = vdarr(*,igood)      ;use good chunks only

; dave - experiment
;for ob=0,24 do begin
;    n_max = n_elements(vdarr[ob,*]) - 1
;    n_min = n_max - 150
;
;    vdarr[ob, n_min:n_max].(veldex) = vdarr[ob, 0:150].(veldex)
;endfor


;
IF not keyword_set(nozero) then begin
    ; Add a CONSTANT to the velocities in each chunk set to     
    ; force the median velocity of each chunk set to be ZERO.
    ; This preserves the velocity variation from observaton to observation.
    ; Designed to account for errors in the values of DSST.w0 for each chunk.
    ; Note: This also removes the zero-pt. calibration of the velocity scale.
;

  if n_elements(sp2) eq 1 then if sp2 eq 1 then veldex = first_el(where(tagnam eq 'SP2'))
  FOR n = 0 , nchunk - 1 do begin         ;Loop thru chunks
     vset = vdarr(*,n).(veldex)           ;rename vel's of chunk set
     vind = where(vdarr(*,n).(fitdex) lt 99.)
     mnchu = mean(vset(vind))             ;mean vel of chunk set
     ;;;replace with median (Fischer 11Feb2013)
     ;mnchu = median(vset(vind))             ;median vel of chunk set
     vdarr(*,n).(veldex) = vset - mnchu   ;shift vels of chunk set
  END

  print,'  VDCLEAN: Set "middle" Velocity of Each Chunk Set to 0'
ENDIF

;
;Within a chunk set, determine the discrepancy between the velocity 
;of each member and the median velocity of its respective observation. 
;
  diff = fltarr(numobs,nchunk)  ;difference: vel - median_obs

  FOR ob = 0,numobs-1 do begin
    medvel = median(vdarr(ob,*).(veldex))       ;median Velocity off obs.
    diff(ob,*) = vdarr(ob,*).(veldex) - medvel  ;Vel_chunk - median_obs
;;;debugging the NaN problem ;;fischer Jan 17,2013
;    print,medvel
;    res=finite(diff[ob,*])  & gg=where(res eq 0,ngg)
;    if ngg gt 0 then print,vdarr[ob].obnm
 ;   xd=where(diff[ob,*] eq 0.0d,nxd)
 ;   if nxd gt 0 then stop
  END

  sigma = fltarr(nchunk)                  ;sigma of vel_chunk - median_obs
;  rv = assign_weights5(vdarr, nchunk, numobs, fitdex, diff, sigma, cut_thres=0.65)  
  rv = assign_weights1(vdarr, nchunk, numobs, fitdex, diff, sigma)  
;  rv=assign_weights2(vdarr, diff,  nchunk, fitdex, numobs) ; bad! 
;  rv=assign_weights3(vdarr,  nchunk, fitdex, numobs)  ;sucks!


;print,minmax(diff)
dum=abs(diff)   &  idum=sort(dum)  &  dum=dum[idum]
max_ind=round(0.97*n_elements(dum))  ;keep up to 97th percentile points
lim_dum=dum[max_ind]  ; 0.99 limit
plothist,diff,xra=[-lim_dum,lim_dum],bin=.5
x=where(abs(diff) gt lim_dum,nx)   ; set outliers to zero weight
if nx gt 0 then vdarr(x).weight=0.
;xx=where(vdarr.fit gt 1.3,nxx)
;if nxx gt 0 then vdarr[xx].weight=0.


 print,'  VDCLEAN: CHUNK WEIGHTS ADJUSTED, according to chunk set scatter.'
;

  ;stop; dave - here do: plot, vdarr[ob,*].weight, psym=4, yr=[0,0.0004]


GOTO, NORMAL

; dave experiment
max_ob_index = 24
array1 = dblarr(100*150)
index = 0
for ob=0,max_ob_index do begin
    for n=0,149 do begin
        array1[index] = vdarr[ob,n].ivel
        index = index + 1
    endfor 
endfor

array2 = dblarr(100*150)
index = 0
for ob=0,max_ob_index do begin
    n_max = n_elements(vdarr[ob,*]) - 1
    n_min = n_max - 149

    for n=n_min,n_max do begin
        array2[index] = vdarr[ob,n].ivel
        index = index + 1
    endfor 
endfor

print, 'dave. array1: mean = ', mean(array1), ', stdev = ', stddev(array1)
print, 'dave. array2: mean = ', mean(array2), ', stdev = ', stddev(array2)

; dave - fix this (thsigma)
;print,''
;print,'    Reject chunk-sets having:'
;print,'           Photon-Limited Error >', fix(therr),' m/s'
;print,'           Median CHISQ >', thchi
;print,'           Scatter from median >',thsigma,' m/s'
print,' Retaining ',strtrim(string(nchunk),2), $
         ' out of ',strtrim(string(orignum),2), ' chunks.'
;print,' '

NORMAL:

;
;IF keyword_set(plot) then begin
;     !p.color=100
;     xtit = '!6 Chunk-Set <Chisq>'
;     plothist,badness,x,y,bin=.1,title='Histogram of <ChiSq>',xtit=xtit
;     oplot,[badness(isort(threshel)),badness(isort(threshel))],[0,100],co=200
;     !p.color=200
;     xyouts,min(x)+0.0*(max(x)-min(x)),0.93*max(y),'GOOD CHIs',size=1.5
;     xyouts,min(x)+0.7*(max(x)-min(x)),0.93*max(y),'BAD GUYs',size=1.5
;     xyouts,min(x)+0.99*(max(x)-min(x)),0.8*max(y),'UGLY',size=0.8
;     !p.color=250
;END
return
end
