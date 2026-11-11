import numpy as np, pickle
from angdss import *
rng=np.random.default_rng(7)
N=256; Ns=32; blk=N*Ns
env=pickle.load(open('env.pkl','rb'))[1]; g,sm2=env['g'],env['sm2']
M20=654729075.0
def qpsk(n):
    b=rng.integers(0,2,(n,2)); return ((1-2*b[:,0])+1j*(1-2*b[:,1]))/np.sqrt(2)
def chips(nsym,long=True,pn=None):
    s=qpsk(nsym)
    c=(2*rng.integers(0,2,(N,nsym))-1) if long else np.repeat(pn[:,None],nsym,1)
    return (c*s[None,:]).T.reshape(-1)
def mask(x,sm):
    return np.real(x)*sm*np.abs(rng.standard_normal(x.size))+1j*np.imag(x)*sm*np.abs(rng.standard_normal(x.size))
def cg(n,v): return np.sqrt(v/2)*(rng.standard_normal(n)+1j*rng.standard_normal(n))
# pooled 20th moment over Kb blocks at -10 dB average SNR
snr=10**(-1.0); nv=1/snr
Kb=2000; idx=np.arange(Kb)  # first 2000 blocks of trace (contiguous)
acc={k:[0.0,0.0,0] for k in ['noise','gdss','angdss','noise_const','gdss_const']}
for k in idx:
    x=chips(Ns)
    w=cg(blk,nv*g[k]); w2=cg(blk,nv*g[k]); wc=cg(blk,nv)
    sig={'noise':w2,'gdss':mask(x,1)+w,'angdss':mask(x,np.sqrt(sm2[k]))+cg(blk,nv*g[k]),'noise_const':wc,'gdss_const':mask(x,1)+cg(blk,nv)}
    for n,y in sig.items():
        r=np.concatenate([y.real,y.imag]); acc[n][0]+=np.sum(r**2); acc[n][1]+=np.sum(r**20); acc[n][2]+=r.size
for n,(s2,s20,cnt) in acc.items():
    m2=s2/cnt; m20=s20/cnt/m2**10
    print('pooled 20th moment dev %-12s %8.1f %%'%(n,abs(m20-M20)/M20*100))
# per-block (single 8192-chip window) deviation, mean over blocks
