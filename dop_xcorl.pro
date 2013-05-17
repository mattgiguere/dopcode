FUNCTION DOP_XCORL, ISS_OBNM, OBSNM, BLUE_DISP

c_light=2.99792458d8

rdsi,iss,iss_obnm
rdsi,obs,obsnm

snip_iss=iss[1000:1300,55]
snip_obs=obs[1000:1300,55]
contf,snip_iss,c_iss,sbin=20,nord=2
contf,snip_obs,c_obs,sbin=20,nord=2

xcorlb,snip_iss/c_iss,snip_obs/c_obs,25,shift
init_z = c_light*shift*blue_disp/xcorl_lambda
