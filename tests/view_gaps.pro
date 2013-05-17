;; This program plots a histogram of the wavelength gaps between chunks in an order.
;; Example:
;;   wav0 = chunk[i].wav0 + chunk[i].wcof[0] + 40*chunk[i].wcof[1]
;;   wav1 = chunk[i+1].wav0 + chunk[i+1].wcof[0]
;;   gap = wav1 - wav0

PRO VIEW_GAPS, vd_filename=vd_filename, vd_tag=vd_tag, vd_dir=vd_dir, vdst=vdst, call_sm_wav=call_sm_wav, call_sm_disp=call_sm_disp, use_simple_slope=use_simple_slope, window_num=window_num, xr=xr

if ~keyword_set(vd_dir) then begin
    vd_dir = '/mir1/files/'
endif

restore,vd_dir+vd_filename

; Lick specific
total_chunks = n_elements(vd)
if total_chunks eq 704 then chunks_per_order = 44
if total_chunks eq 352 then chunks_per_order = 22
num_orders = total_chunks / chunks_per_order;   16 for Lick
start_order = 38

; comment out specially for 'ra' observations
;total_chunks = 272
;chunks_per_order = 17
;num_orders = total_chunks / chunks_per_order;   16 for Lick
;start_order = 6


if keyword_set(vd_filename) then begin
    vd_names = strarr(1)
    vd_names[0] = vd_dir + vd_filename
endif else begin
   if (keyword_set(vd_tag)) then begin
    search = vd_dir + vd_tag + '*'
    vd_names = findfile(search)
   endif 
endelse
if (keyword_set(vdst)) then begin
    vd = vdst
endif

if (~keyword_set(vdst)) then begin
    gaps = dblarr(n_elements(vd_names) * (total_chunks - 1))
endif else begin
    vd_names = strarr(1)
    gaps = dblarr(1 * (total_chunks - 1))
endelse
gap_index = long(0)

; Go through vd files
for v=0,n_elements(vd_names)-1 do begin

    if (~keyword_set(vdst)) then begin
        print, 'Measuring gaps for vd file: ', vd_names[v]
        restore, vd_names[v]
    endif

    if (keyword_set(call_sm_disp)) then begin
        print, 'calling sm_disp'
        sm_disp, vd, nvd
        vd = nvd
    endif
    

    if (keyword_set(call_sm_wav)) then begin
        print, 'calling sm_wav'
        if (keyword_set(use_simple_slope)) then begin
            sm_wav, vd, wc, pord=4, /simple_slope
        endif else begin
            sm_wav, vd, wc, pord=4
        endelse
    endif

    ; Go through the orders.
    for ordr=0,num_orders-1 do begin
        start_chunk = ordr * chunks_per_order
        end_chunk = start_chunk + chunks_per_order - 1

        ; Go through all of the chunks in this order, starting with the second
        for i=start_chunk+1,end_chunk do begin
            prev_chunk = vd[i-1]
            this_chunk = vd[i]
            num_pixels = this_chunk.pixob - prev_chunk.pixob
            wav0 = double(prev_chunk.w0) + double(prev_chunk.wcof[0]) + (num_pixels)*double(prev_chunk.wcof[1])
            wav1 = double(this_chunk.w0) + double(this_chunk.wcof[0])
            gaps[gap_index] = wav1 - wav0
            gap_index = gap_index + 1
        endfor; chunks
    endfor ; orders

endfor ;  vds

max_gap = max(gaps)
min_gap = min(gaps)
mean_gap = mean(gaps)
stddev_gap = stddev(gaps)

print, 'Calculated gaps for ', n_elements(gaps), ' chunk pairs.'
print, 'largest gap: ', max_gap
print, 'smallest gap: ', min_gap
print, 'mean gap: ', mean_gap
print, 'stddev: ', stddev_gap

if ~keyword_set(xr) then begin
    xval = max([abs(min_gap), abs(max_gap)])
    xrange = [-xval, xval]
endif else begin
    xrange = xr
endelse

; White background works better for printing. - dave 2/24/06
nicecolor, nc
!p.background = nc.white
!p.color = nc.black

xrange=[-0.0002,0.0002]

title='Histogram of chunk gaps for' + strcompress(string(n_elements(gaps))) + ' chunk pairs'

plothist, gaps, bin=0.000001, title=title, xrange=xrange, xtitle='Gap between neighboring chunks (Ansgtroms)', ytitle='Number of chunks'

; If a specifc vd file/structure was given, plot the gaps for each order.
if (keyword_set(vd_filename) or keyword_set(vdst)) then begin
    if ~keyword_set(window_num) then begin
        window_num = 1
    endif
    window, window_num, XSize=800, YSize=800

    !P.Multi = [0, 4, 4, 0, 0]

    ; Go through the orders.
    for ordr=0,num_orders-1 do begin
        start_chunk = (ordr * (chunks_per_order-1))
        end_chunk = start_chunk + (chunks_per_order-1) - 1
        x = dindgen(end_chunk - start_chunk + 1)
        x = x + start_chunk 
        ytitle = 'wvln gap'
        xtitle = 'order' + strcompress(string(ordr+start_order))
        ;yrange = max([abs(min(gaps)), abs(max(gaps))])
        plot, x, gaps[start_chunk:end_chunk - 1], xtitle=xtitle, ytitle=ytitle
    endfor

endif

END
