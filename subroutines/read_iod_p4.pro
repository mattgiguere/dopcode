pro read_iod_p4, cell, w,s, w1, w2, pad=pad

dir='/home/fischer/dop2/atlas/'
if ~keyword_set(cell) then cell = 'iodine_p4_pnnl_lo_apod.dat'
file = dir+cell
if ~keyword_set(pad) then pad = 1.

restore, file ; wvac, tran
x=where(wvac gt w1-pad and wvac lt w2+pad,nx)
w=wvac[x]
s=tran[x]

end ;pro
