import numpy as np, time, pickle
from angdss import *
t=time.time()
N=256
tbl=Table(N)
print('table',time.time()-t)
sa=np.arange(-30,5.0001,0.1)
print('static GDSS 1e-5:',required_snr(sa,tbl(10**(sa/10))))
res={}
for seed in [1,2,3]:
    rng=np.random.default_rng(seed)
    g=noise_environment(rng,100000,3)
    s2=estimate_noise_var(rng,g,256,'median',0.01,20)
    sm2=mask_power(s2); smg=mask_power(g)
    C=ber_curves(sa,g,sm2,N,tbl); Cg=ber_curves(sa,g,smg,N,tbl)
    r={k:required_snr(sa,C[k]) for k in C}; r['genie']=required_snr(sa,Cg['angdss'])
    print(seed,{k:round(v,2) for k,v in r.items()}, 'gain',round(r['gdss']-r['angdss'],2), 'snrstdA',np.std(10*np.log10(sm2/g)).round(2),'esterr',np.std(10*np.log10(s2/g)).round(2))
    res[seed]=dict(g=g,s2=s2,sm2=sm2,smg=smg)
pickle.dump((tbl.snr_db,tbl.lb),open('tbl256.pkl','wb'))
pickle.dump(res,open('env.pkl','wb'))
