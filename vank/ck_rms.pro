pro ck_rms

restore,'vstbank/vsthip8102.dat'

narr=25
rms=fltarr(narr)  


n1=where(cf1.jd gt 14373. and cf1.jd lt 14374.,nn1)  & d1=cf1(n1) &  rms(0)=stdev(d1.mnvel)
n2=where(cf1.jd gt 14374. and cf1.jd lt 14375.,nn2)  & d2=cf1(n2) &  rms(1)=stdev(d2.mnvel)
n3=where(cf1.jd gt 14381. and cf1.jd lt 14382.,nn3)  & d3=cf1(n3) &  rms(2)=stdev(d3.mnvel)
n4=where(cf1.jd gt 14382. and cf1.jd lt 14383.,nn4)  & d4=cf1(n4) &  rms(3)=stdev(d4.mnvel)
n5=where(cf1.jd gt 14395. and cf1.jd lt 14396.,nn5)  & d5=cf1(n5) &  rms(4)=stdev(d5.mnvel)
n6=where(cf1.jd gt 14396. and cf1.jd lt 14397.,nn6)  & d6=cf1(n6) &  rms(5)=stdev(d6.mnvel)
n7=where(cf1.jd gt 14397. and cf1.jd lt 14398.,nn7)  & d7=cf1(n7) &  rms(6)=stdev(d7.mnvel)
n8=where(cf1.jd gt 14398. and cf1.jd lt 14399.,nn8)  & d8=cf1(n8) &  rms(7)=stdev(d8.mnvel)
n9=where(cf1.jd gt 14399. and cf1.jd lt 14400.,nn9)  & d9=cf1(n9) &  rms(8)=stdev(d9.mnvel)
n10=where(cf1.jd gt 14407. and cf1.jd lt 14408.,nn10)  & d10=cf1(n10) &  rms(9)=stdev(d10.mnvel)
n11=where(cf1.jd gt 14408. and cf1.jd lt 14409.,nn11)  & d11=cf1(n11) &  rms(10)=stdev(d11.mnvel)
n12=where(cf1.jd gt 14409. and cf1.jd lt 14410.,nn12)  & d12=cf1(n12) &  rms(11)=stdev(d12.mnvel)
n13=where(cf1.jd gt 14414. and cf1.jd lt 14415.,nn13)  & d13=cf1(n13) &  rms(12)=stdev(d13.mnvel)
n14=where(cf1.jd gt 14416. and cf1.jd lt 14417.,nn14)  & d14=cf1(n14) &  rms(13)=stdev(d14.mnvel)
n15=where(cf1.jd gt 14424. and cf1.jd lt 14425.,nn15)  & d15=cf1(n15) &  rms(14)=stdev(d15.mnvel)
n16=where(cf1.jd gt 14425. and cf1.jd lt 14426.,nn16)  & d16=cf1(n16) &  rms(15)=stdev(d16.mnvel)
n17=where(cf1.jd gt 14426. and cf1.jd lt 14427.,nn17)  & d17=cf1(n17) &  rms(16)=stdev(d17.mnvel)
n18=where(cf1.jd gt 14427. and cf1.jd lt 14428.,nn18)  & d18=cf1(n18) &  rms(17)=stdev(d18.mnvel)
n19=where(cf1.jd gt 14428. and cf1.jd lt 14429.,nn19)  & d19=cf1(n19) &  rms(18)=stdev(d19.mnvel)
n20=where(cf1.jd gt 14433. and cf1.jd lt 14434.,nn20)  & d20=cf1(n20) &  rms(19)=stdev(d20.mnvel)
n21=where(cf1.jd gt 14434. and cf1.jd lt 14435.,nn21)  & d21=cf1(n21) &  rms(20)=stdev(d21.mnvel)
n22=where(cf1.jd gt 14435. and cf1.jd lt 14436.,nn22)  & d22=cf1(n22) &  rms(21)=stdev(d22.mnvel)
n23=where(cf1.jd gt 14447. and cf1.jd lt 14448.,nn23)  & d23=cf1(n23) &  rms(22)=stdev(d23.mnvel)
n24=where(cf1.jd gt 14448. and cf1.jd lt 14449.,nn24)  & d24=cf1(n24) &  rms(23)=stdev(d24.mnvel)
n25=where(cf1.jd gt 14449. and cf1.jd lt 14450.,nn25)  & d25=cf1(n25) &  rms(24)=stdev(d25.mnvel)


for i=0,narr-1 do print,i+1, rms(i)
print,mean(rms)


stop
end
