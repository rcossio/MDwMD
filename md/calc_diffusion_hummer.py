from Dfit import Dcov


res = Dcov(m=100,fz='tmp.xvg',tmin=1,tmax=20,nseg=20, dt=10, nitmax=50000)
res.run_Dfit()
res.analysis(tc=50)

#import Dfit

#res = Dfit.Dfit.Dcov(m=99,fz='tmp.xvg',tmin=1,tmax=10,nseg=20, dt=10, nitmax=300)
#res.run_Dfit()
#res.analysis(tc=20)