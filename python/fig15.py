import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
o=json.load(open('extra.json'))
B=(0,0.45,0.74); O=(0.85,0.33,0.10)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.linewidth':0.8})
fig,ax=plt.subplots(1,2,figsize=(12.5,5.4))
D=np.array(o['D']); c=o['curves']
a=ax[0]
a.plot(D,c['GDSS_1024']['q95'],'-',color=B,lw=2,label='GDSS, $n_w$ = 1024')
a.plot(D,c['ANGDSS_1024']['q95'],'-',color=O,lw=2.2,label='ANGDSS, $n_w$ = 1024')
a.plot(D,c['GDSS_8192']['q95'],'--',color=B,lw=1.6,label='GDSS, $n_w$ = 8192 (full block)')
a.plot(D,c['ANGDSS_8192']['q95'],'--',color=O,lw=1.8,label='ANGDSS, $n_w$ = 8192 (full block)')
a.axhline(0.2,color='r',ls=':',lw=1.3)
a.text(0.4,0.215,'design limit $P_D$ = 0.2',color='r',fontsize=9)
a.set_xlim(0,25);a.set_ylim(0,1.02);a.grid(alpha=.35)
a.set_xlabel('Warden SNR offset $\\Delta$ below the receiver (dB)');a.set_ylabel('95th-percentile $P_D$ over blocks')
a.set_title('(a) Detectability at equal reliability (BER = $10^{-5}$)',fontweight='bold',fontsize=11.5)
a.legend(loc='upper center',bbox_to_anchor=(0.5,-0.17),ncol=2,fontsize=9,framealpha=1)
b=ax[1]
taus=sorted(int(k) for k in o['causal']); req=[o['causal'][str(t)]['req'] for t in taus]
b.semilogx(taus,req,'-d',color=O,lw=2,mfc=O,label='ANGDSS, causal controller')
b.axhline(o['offline'],color=O,ls=':',lw=1.6,label='ANGDSS, block-average controller')
b.axhline(o['op']['GDSS'],color=B,lw=1.6,label='Fixed-mask GDSS')
b.set_xlim(15,2500);b.set_ylim(-9.6,-4.4);b.grid(alpha=.35,which='both')
b.set_xlabel('Averaging time constant $\\tau$ (blocks)');b.set_ylabel('Average SNR for BER = $10^{-5}$ (dB)')
b.set_title('(b) Causal mask controller, $\\sigma$ = 3 dB',fontweight='bold',fontsize=11.5)
b.legend(loc='upper center',bbox_to_anchor=(0.5,-0.17),ncol=2,fontsize=9,framealpha=1)
for s in ax:
    for sp in s.spines.values(): sp.set_visible(True)
plt.tight_layout(); plt.savefig('results/fig15_equal_ber_causal.png',dpi=170)
