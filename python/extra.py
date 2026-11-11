import numpy as np, pickle, json
from angdss import *
N=256
tbl=Table(N); sa=np.arange(-30,5.0001,0.1)
# A: realisation spread
gains=[];reqA=[];reqG=[];reqD=[]
for seed in range(1,11):
    rng=np.random.default_rng(100+seed)
    g=noise_environment(rng,100000,3); s2=estimate_noise_var(rng,g,256,'median',0.01,20); sm2=mask_power(s2)
    C=ber_curves(sa,g,sm2,N,tbl)
    a,gg,d=required_snr(sa,C['angdss']),required_snr(sa,C['gdss']),required_snr(sa,C['dsss'])
    reqA.append(a);reqG.append(gg);reqD.append(d);gains.append(gg-a)
out={'spread':dict(gain_mean=np.mean(gains),gain_std=np.std(gains,ddof=1),gain_min=min(gains),gain_max=max(gains),
     A_mean=np.mean(reqA),A_std=np.std(reqA,ddof=1),G_mean=np.mean(reqG),G_std=np.std(reqG,ddof=1),D_mean=np.mean(reqD),D_std=np.std(reqD,ddof=1))}
print(out['spread'])
# use seed-1 env (matches main numbers)
env=pickle.load(open('env.pkl','rb'))[1]; g,s2,sm2=env['g'],env['s2'],env['sm2']
C=ber_curves(sa,g,sm2,N,tbl)
op={'GDSS':required_snr(sa,C['gdss']),'ANGDSS':required_snr(sa,C['angdss'])}
pw={'GDSS':np.ones_like(g),'ANGDSS':sm2}
D=np.arange(0,25.01,0.25)
curves={}
for nw in [1024,8192]:
    for s in ['GDSS','ANGDSS']:
        q95=[];mean=[]
        for d in D:
            gam=10**((op[s]-d)/10)*pw[s]/g
            p=pd_radiometer(gam,nw)
            q95.append(np.quantile(p,0.95)); mean.append(p.mean())
        q95=np.array(q95);mean=np.array(mean)
        dreq=float(np.interp(0.2,q95[::-1],D[::-1]))
        curves[f'{s}_{nw}']=dict(q95=q95.tolist(),mean=mean.tolist(),dreq=dreq,pd0_mean=float(mean[0]),pd0_q95=float(q95[0]))
        print(s,nw,'op',round(op[s],2),'PD at D=0 mean %.3f q95 %.3f'%(mean[0],q95[0]),'Delta for q95<=0.2: %.2f dB'%dreq)
# KL bound for ANGDSS at its Delta_req, nw=1024 (median block)
for nw in [1024,8192]:
    d=curves[f'ANGDSS_{nw}']['dreq']; gam=10**((op['ANGDSS']-d)/10)*sm2/g
    Dkl=nw*(gam-np.log1p(gam)); lb=1-np.sqrt(np.median(Dkl)/2)
    print('KL bound nw',nw,'median D',np.median(Dkl).round(4),'P_FA+P_MD >=',lb.round(3))
    curves[f'ANGDSS_{nw}']['kl_lb']=float(lb)
# B: causal controller
caus={}
for tau in [20,50,100,200,500,1000,2000]:
    smc=mask_power_causal(s2,tau=tau)
    Cc=ber_curves(sa,g,smc,N,tbl)
    r=required_snr(sa,Cc['angdss'])+10*np.log10(smc.mean())
    caus[tau]=dict(req=float(r),meanpow_dB=float(10*np.log10(smc.mean())),clipfrac=float(np.mean(smc>=4-1e-9)))
    print('tau',tau,caus[tau])
out.update(op=op,D=D.tolist(),curves=curves,causal=caus,offline=float(op['ANGDSS']))
json.dump(out,open('extra.json','w'),indent=1)
