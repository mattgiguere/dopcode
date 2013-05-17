FUNCTION DOP_VDARR_SETUP, OBSNM, DOPENV, NCHUNKS

npix=dopenv.n_pix
npar_psf=n_elements(dopenv.psfpix)

vd={obnm: obsnm,     $ ;observation name
        ordt: 0,         $ ;order for ISS
        ordob:0,         $ ;order of observation
        pixt:0,          $ ;starting pixel for chunk of ISS
        pixob:0,         $ ;starting pixel for observation
        w0:0,            $ ;
        wcof:fltarr(2),  $ ;wavelength coef's for chunk
        cts:0L,          $ ;photon counts 
        z:0.0d,          $ ;doppler z (no BC correction)
        fit:0.0,         $ ;chisq fit
        npix:npix,       $ ;number of pixels in chunk
        vel: 0.0d,        $ ;z*c_light 
        weight:0.0,      $ ;chunk weight
        psf:fltarr(121), $   ;amplitudes of psf
        par:fltarr(npar_psf)}
vd=replicate(vd,nchunks)

return,vd

end
