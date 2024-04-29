from Dfit import Dcov


#res = Dcov(m=100,fz='tmp.xvg',tmin=1,tmax=20,nseg=20, dt=10, nitmax=50000) #5PTI
#res = Dcov(m=79,fz='tmp.xvg',tmin=1,tmax=25,nseg=25, dt=10, nitmax=50000) #6LYZ
#res = Dcov(m=24,fz='tmp.xvg',tmin=1,tmax=100,nseg=20, dt=10, nitmax=50000) #2KV4
res = Dcov(m=24,fz='tmp.xvg',tmin=1,tmax=100,nseg=20, dt=10, nitmax=50000) #2KV4b
#res = Dcov(m=20,fz='tmp.xvg',tmin=1,tmax=100,nseg=20, dt=10, nitmax=50000) #4F5S


res.run_Dfit()

#res.analysis(tc=50) #5PTI
#res.analysis(tc=130) #6LYZ
#res.analysis(tc=400) #2KV4
res.analysis(tc=200) #2KV4b


#import Dfit

#res = Dfit.Dfit.Dcov(m=99,fz='tmp.xvg',tmin=1,tmax=10,nseg=20, dt=10, nitmax=300)
#res.run_Dfit()
#res.analysis(tc=20)
