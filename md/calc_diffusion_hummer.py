from Dfit import Dcov


#res = Dcov(m=100,fz='tmp.xvg',tmin=1,tmax=20,nseg=20, dt=10, nitmax=50000) #5PTI
res = Dcov(m=20,fz='tmp.xvg',tmin=1,tmax=12,nseg=20, dt=10, nitmax=50000) 
#res = Dcov(m=79,fz='tmp.xvg',tmin=1,tmax=25,nseg=25, dt=10, nitmax=50000) #6LYZ

res.run_Dfit()

res.analysis(tc=50) #5PTI
#res.analysis(tc=130) #6LYZ

#import Dfit

#res = Dfit.Dfit.Dcov(m=99,fz='tmp.xvg',tmin=1,tmax=10,nseg=20, dt=10, nitmax=300)
#res.run_Dfit()
#res.analysis(tc=20)
