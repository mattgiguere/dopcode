pro dr_psf_smooth, order=order, pix=pix

if ~keyword_set(order) then order = 42
if ~keyword_set(pix) then pix = 900

restore,'/mir1/files/vdaHIP8102_rg92.234'
psfav_file='vdaHIP8102_rg92.234'
chunk=chunk_arr
dopenv=dop_init('rg92.234', 'dsstHIP8102q_rg92.dat', 4494.078)
ip_av=dop_psf_smooth(chunk, order, pix, dopenv=dopenv, psfav_file = psfav_file, /plot) 

stop
end
