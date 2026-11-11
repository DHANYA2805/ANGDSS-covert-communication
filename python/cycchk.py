import numpy as np, pickle
from angdss import *
rng=np.random.default_rng(11)
N=256; Ns=32; blk=N*Ns
env=pickle.load(open('env.pkl','rb'))[1]; g,sm2=env['g'],env['sm2']
def ees(x,Nw):
    K=x.size//Nw; X=np.fft.fft(x[:K*Nw].reshape(K,Nw).T,axis=0); P=np.mean(np.abs(X)**2,1)
    M=Nw//2-1; out=np.zeros(M)
    for m in range(1,M+1):
        Xs=np.roll(X,-m,0); Ps=np.roll(P,-m); Sa=np.mean(Xs*np.conj(X),1); out[m-1]=np.mean(np.abs(Sa)/np.sqrt(Ps*P))
    return np.arange(1,M+1)/Nw,out,np.sqrt(np.pi/(4*K))
nsym=4096; Kb=nsym//Ns
s=((1-2*rng.integers(0,2,nsym))+1j*(1-2*rng.integers(0,2,nsym)))/np.sqrt(2)
pnS=2*rng.integers(0,2,N)-1; pnL=2*rng.integers(0,2,(N,nsym))-1
cS=(pnS[:,None]*s[None,:]).T.reshape(-1); cL=(pnL*s[None,:]).T.reshape(-1)
perchip=lambda v: np.repeat(v[:Kb],blk)
def mask(x,sm): return np.real(x)*sm*np.abs(rng.standard_normal(x.size))+1j*np.imag(x)*sm*np.abs(rng.standard_normal(x.size))
snr=10**(-1.2)
gC=perchip(g); aC=perchip(sm2)
w=lambda: np.sqrt(gC/2)*(rng.standard_normal(gC.size)+1j*rng.standard_normal(gC.size))
Nw=4*N
for name,x in [('noise only (drift)',0*cL),('ANGDSS long code',mask(cL,np.sqrt(aC))),('ANGDSS repeating code',mask(cS,np.sqrt(aC)))]:
    y=np.sqrt(snr)*x+w()
    a,e,fl=ees(y,Nw); m=np.round(a*Nw).astype(int); on=(m%(Nw//N)==0)&(m>Nw//N)
    print('%-24s excess at k/N: %.3f   max off-line excess %.3f'%(name,e[on].mean()/fl-1,(e[~on]/fl-1).max()))
