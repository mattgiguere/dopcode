pro view_psf,vd_array=vd_array,vd_filename=vd_filename,original=original,contour=contour,del_pix=del_pix,del_ord=del_ord,start_window=start_window,best_fit=best_fit,median_fit=median_fit,mf_id=mf_id,make_jpegs=make_jpegs,psfpix=psfpix,psfsig=psfsig

; Declare common block needed for gpfunc.pro.
;common psfstuff,sig,pix,obpix
;common psfstuff,psfsig,psfpix,obpix,param

osamp=4
fwhm_osamp = 10
hm_delta = 0.005
orddist=15; for Lick

if not keyword_set(start_window) then begin
    start_window=0
endif
next_window = start_window

; Initialize the set of Gaussians used to model the PSF.
sigpsf=fltarr(11)
pixpsf=fltarr(11)

if (keyword_set(psfpix) and keyword_set(psfsig)) then begin
    pixpsf = psfpix
    sigpsf = psfsig

endif else begin
    print, 'NOTICE: Using PSF description hard-coded in file! Make SURE this is the same description used by crank to generated your vd(s)!'
    wait, 1

    ; PSF RF L
;    pixpsf=1.1*[0.00,-3.25,-2.60,-1.95,-1.30, -0.65, 0.65, 1.30, 1.95, 2.60, 3.25, 0.0, 0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
;    sigpsf=-0.2+[0.9, 0.40, 0.40, 0.40, 0.40,  0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.0, 0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ;pixpsf=1.3*[0.00,-3.25,-2.60,-1.95,-1.30, -0.65, 0.65, 1.30, 1.95, 2.60, 3.25, 0.0, 0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ;sigpsf=[1.0, 0.40, 0.40, 0.40, 0.40,  0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.0, 0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ; PSF description cooked up by Dr. Fischer and I on 3/10/06 for cranking lots of data for 9826 (crossing multiple CCDs).
;    pixpsf=[0.00,-2.50,-1.80, -1.10,-0.80, -0.50, 0.50, 0.80, 1.10, 1.80, 2.50, 0.0, 0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
;    sigpsf=[0.9, 0.40, 0.40, 0.40, 0.40,  0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.0, 0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ; PSF used by Pauls' code at UCB
     pixpsf=[0.00,-2.40,-2.10,-1.60,-1.10,-0.60, 0.60, 1.10, 1.60, 2.10, 2.40, 0.00, 0.00, 0.00,0.00]
     sigpsf=[0.40, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.00, 0.00, 0.00,0.00]

;jj's subaru psf
;        psfpix=[0.0, -2.0, 2.0, -4.0, 4.0]
;        psfsig=[2.3,  1.5, 1.5,  1.5, 1.5]

; PSF description ************ tag = f
;    pix = [0.00,-2.70,-2.10,-1.50,-1.00,-0.50, 0.50, 1.00, 1.50, 2.10, 2.70, 0.00, 0.00, 0.00, 0.00]
;    sig = [0.85, 0.60, 0.60, 0.60, 0.60, 0.40, 0.40, 0.60, 0.60, 0.60, 0.60, 0.00, 0.00, 0.00, 0.00]

; PSF description ************ tag = g
;    pix = [0.00,-1.80,-1.30,-0.90,-0.60,-0.30, 0.30, 0.60, 0.90, 1.30, 1.80, 0.00, 0.00, 0.00, 0.00]
;    sig = [0.9,  0.60, 0.60, 0.50, 0.40, 0.40, 0.40, 0.40, 0.50, 0.60, 0.60, 0.00, 0.00, 0.00, 0.00]



; Actually, this is for JJ EMU stuff. Comes from his rdbcvel, which is called
; from his cf.pro, since I'm not presently specifying a psf description in
; my dop_driver.pro code for lick.
;pixpsf = [0.00000, -3.70000, -3.00000, -2.30000, -1.60000, -0.900000, 0.900000, 1.60000, 2.30000, 3.00000, 3.70000, 0.00000, 0.00000, 0.00000, 0.00000]
;sigpsf = [0.800000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000,  0.00000,  0.00000,  0.00000,  0.00000]


endelse

sig=sigpsf; For gpfunc.
pix=pixpsf; For gpfunc.

if keyword_set(vd_filename) then begin
    restore, vd_filename
    vd0 = vd
    vd_array = 0
endif else begin
    if keyword_set(vd_array) then begin
        vd0 = vd_array[0:703]
    endif else begin
        print, 'Either vd_filename of vd_array must be set'
        return
    endelse
endelse

; Create a 2D array that will be used for plotting
num_orders = n_elements(unique(vd0.ordob))
num_chunks_per_order =  n_elements(vd0) / num_orders
ccd_chunks = fltarr(num_chunks_per_order, num_orders)
order_axis = indgen(num_orders)
order_axis = order_axis + min(vd0.ordob)
chunk_axis = indgen(num_chunks_per_order)

; Set del_pix and del_ord for psfav
if not keyword_set(del_ord) then del_ord=3
if not keyword_set(del_pix) then del_pix=60

if not keyword_set(original) then begin
    print, 'psfav will be called with del_ord=', strcompress(string(del_ord)), $
           ', and del_pix=', strcompress(string(del_pix))
end

; For each chunk of the CCD, obtain either:
;  - The PSF calculated for that chunk (could be from the first or second call
;    to crank, depending on the input file). Requires the 'original' key.
;  - The average PSF for that chunk, given by the call to psfav and the
;    parameters del_ord, del_pix. This is the default operation.
;
; Once that's been determined, calculate the FWHM for the PSF.

numpix = 15  ; PSF is modeled b/n -numpix <= 0 <= +numpix pixels.
numpoints = ((2 * numpix) * osamp) + 1; Oversampling by factor of osamp. 

; Create the oversampled x axis (-numpix <= 0 <= +numpix).
x = findgen(numpoints)
x = x/osamp  ; set scale of range
x = x - numpix  ; change range to span from -numpix to +numpix.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; best_fit code ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (keyword_set(best_fit)) then begin
    if (vd_array eq 0) then begin
        vd1 = vd0
    endif else begin
        vd1 = vd_array
    endelse

    fit = MIN(vd1.fit, min_ind)
    psf = gpfunc(x, vd1[min_ind].iparam)
    title = 'Original PSF with best fit: chunk ' + strcompress(string(min_ind)) + ', fit =' + strcompress(string(vd1[min_ind].fit)) + ' (no averaging)'
    window, next_window

    plot, x, psf, ps=8, title=title, xr=[-4,4]

    if (keyword_set(make_jpegs)) then begin
        tv = tvrd()
        write_jpeg, 'psf_best_fit_noavg.jpg', tv
    endif 

    next_window = next_window + 1

    ifit = MIN(vd1.ifit, min_ind)
    if (ifit ne 0) then begin ; Check for vdiod files. Only one call to stargrind, so no ifit.
        chunk = vd1[min_ind]
        psfav,vd1,chunk.ordob,chunk.pixob,osamp,sigpsf=sigpsf,pixpsf=pixpsf,psf,del_ord=del_ord,del_pix=del_pix,orddist=orddist
        title = 'Averaged PSF with best fit: chunk ' + strcompress(string(min_ind)) + ', fit =' + strcompress(string(chunk.ifit))
        window, next_window

        plot, x, psf, ps=8, title=title, xr=[-4,4]
        if (keyword_set(make_jpegs)) then begin
            tv = tvrd()
            write_jpeg, 'psf_best_fit_withavg.jpg', tv
        endif 

        next_window = next_window + 1
    endif 
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; median_fit code ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (keyword_set(median_fit)) then begin
    if (vd_array eq 0) then begin
        vd1 = vd0
    endif else begin
        vd1 = vd_array
    endelse

    fit = median(vd1.fit)
    median_ind_array = where(vd1.fit eq fit) 
    median_ind = median_ind_array[0]

    psf = gpfunc(x, vd1[median_ind].iparam)
    title = 'Original PSF with median fit: chunk ' + strcompress(string(median_ind)) + ', fit =' + strcompress(string(vd1[median_ind].fit)) + ' (no averaging)'
    window, next_window

    plot, x, psf, ps=8, title=title, xr=[-4,4]

    if (keyword_set(make_jpegs)) then begin
        tv = tvrd()
        write_jpeg, 'psf_median_fit_noavg.jpg', tv
    endif 


    next_window = next_window + 1

    ifit = median(vd1.ifit)
    median_ind_array = where(vd1.ifit eq ifit) 
    median_ind = median_ind_array[0]

    if (keyword_set(mf_id)) then begin
        mf_id = median_ind
    endif 

    if (ifit ne 0) then begin ; Check for vdiod files. Only one call to stargrind, so no ifit.
        chunk = vd1[median_ind]
        psfav,vd1,chunk.ordob,chunk.pixob,osamp,sigpsf=sigpsf,pixpsf=pixpsf,psf,del_ord=del_ord,del_pix=del_pix,orddist=orddist
        title = 'Averaged PSF with median fit: chunk ' + strcompress(string(median_ind)) + ', fit =' + strcompress(string(chunk.ifit))
        window, next_window

        plot, x, psf, ps=8, title=title, xr=[-4,4]
        if (keyword_set(make_jpegs)) then begin
            tv = tvrd()
            write_jpeg, 'psf_median_fit_withavg.jpg', tv
        endif 
 
        next_window = next_window + 1
    endif 
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



for chunk_index=0,n_elements(vd0)-1 do begin
    chunk = vd0[chunk_index]

    ; Obtain the PSF for this chunk.
    if (not keyword_set(original)) then begin
        ; Calculate a PSF for this chunk which is an average of it's own PSF, and neighboring chunks. This
        ; is what crank does on it's second call, when the /avpsf keyword is set.

        if (keyword_set(vd_filename)) then begin
            ; make a vd array of one element
            vd_array = [vd0]
        endif

        psfav,vd_array,chunk.ordob,chunk.pixob,osamp,sigpsf=sigpsf,pixpsf=pixpsf,psf,del_ord=del_ord,del_pix=del_pix,orddist=orddist

    endif else begin
        ; Original PSF. Don't average, just look at the PSF currently in the vd file.
        psf = gpfunc(x, chunk.iparam)

    endelse

    ; Now that we have this chunk's PSF, calculate the FWHM of it. To do this well, we'll need to 
    ; "fill in the blanks" a bit with fspline, and also more closely spaced points on the x-axis.
    ; (Note to self: I would not need fwhm_osamp if osamp in the program were larger. Might want to
    ; consider making this change, as only about 40-50 of the 121 PSF points are non-zero).
    fwhm_x = findgen(numpoints * fwhm_osamp)
    fwhm_x = fwhm_x / ((numpoints * fwhm_osamp) / (2. * numpix))
    fwhm_x = fwhm_x - numpix

    fwhm_psf = fspline(x, psf, fwhm_x)

    ; Plot the PSF of each chunk as we wizz through the CCD
    ; if keyword_set(original) then begin
    ;     title = 'PSF for chunk' + strcompress(string(chunk_index)) + ', fit=' + strcompress(string(chunk.fit))
    ; endif else begin
    ;        title = 'Averaged PSF for chunk' + strcompress(string(chunk_index)) + ', ifit=' + strcompress(string(chunk.ifit))
    ; endelse
    ; if (chunk_index eq 0) then window, next_window
    ; plot, fwhm_x, fwhm_psf, ps=8, title=title, xr=[-4,4]

    ; Find the Half-Max, and the location of it
    hm = max(fwhm_psf, center)/2.

    ; Find the x values where psf(x) is closest to hm
    delta = hm_delta
    repeat begin ; find value to the right of x=0
        hm_points = where(fwhm_psf le hm+delta and fwhm_psf ge hm-delta)
        rt_points = lonarr(1)
        rt_points[0] = -1
        if hm_points[0] ne -1 then begin
            rt_points = where(hm_points ge center)
        endif
        if hm_points[0] ne -1 and rt_points[0] ne -1 then begin
           rt_points = fwhm_x[hm_points[rt_points]]
           x_right = median(rt_points)
        endif else begin
           delta = delta + hm_delta
        endelse
    endrep until (rt_points[0] ne -1)

    delta = hm_delta
    repeat begin ; find value to the left of x=0
        hm_points = where(fwhm_psf le hm+delta and fwhm_psf ge hm-delta)
        lt_points = lonarr(1)
        lt_points[0] = -1
        if hm_points[0] ne -1 then begin
            lt_points = where(hm_points le center)
        endif
        if hm_points[0] ne -1 and lt_points[0] ne -1 then begin
           lt_points = fwhm_x[hm_points[lt_points]]
           x_left = median(lt_points)
        endif else begin
           delta = delta + hm_delta
        endelse
    endrep until (lt_points[0] ne -1)

    fwhm = x_right - x_left

    ;print, 'fwhm = ', fwhm, ' for chunk ', chunk_index 
    ccd_chunks[(chunk_index mod num_chunks_per_order), (chunk_index / num_chunks_per_order)] = fwhm

endfor; each chunk

;TvLCT, [70,255,0], [70,255,255], [70,0,0], 1
next_window = next_window + 1

window,next_window

ptitle='FWHM of PSFs'
if not keyword_set(original) then begin
    ptitle = ptitle + ', del_ord=' + strcompress(string(del_ord)) + ', del_pix=' + strcompress(string(del_pix))
endif else begin
    ptitle = ptitle + ', original PSFs (no averaging)'
endelse

shade_surf, ccd_chunks, chunk_axis, order_axis, xtitle='Chunk', ytitle='Order', ztitle='PSF FWHM', title=ptitle, zrange=[1, 3]

if (keyword_set(make_jpegs)) then begin
    tv = tvrd()
    write_jpeg, 'psf_fwhm_over_ccd.jpg', tv
endif

if keyword_set(contour) then begin
    next_window = next_window + 1
    window,next_window
    numlevels=12
    contour, ccd_chunks, chunk_axis, order_axis, xtitle='Chunk', ytitle='Order', ztitle='PSF FWHM', $
             title=ptitle, xstyle=1, ystyle=1, nlevels=numlevels, /follow, /fill
    ; Add contour lines on top of the filled conttour shapes.
    contour, ccd_chunks, chunk_axis, order_axis, xtitle='Chunk', ytitle='Order', ztitle='PSF FWHM', $
             title=ptitle, xstyle=1, ystyle=1, nlevels=numlevels, /follow, /overplot
endif

end



