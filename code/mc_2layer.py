#!/usr/bin/env python3
"""
Henyey-Greenstein (anisotropic, g=0.9) white Monte Carlo for the TWO-LAYER slab,
to compute an accurate cortical sensitivity fraction f_cortex = <L_cortex>/<L_total>
for the paper's two-layer model.  (mc_csf.py used g=0; this uses realistic forward
scattering, mu_s = mu_s'/(1-g), which gives accurate absolute partial pathlengths.)
"""
import numpy as np, json, time

def run(geom, N, seed, sds_centers, g=0.9, half_width=2.5,
        z_max=80.0, L_max=500.0, max_iter=20000, n_tissue=1.4):
    rng = np.random.default_rng(seed)
    bounds = np.array(geom['bounds'], float)
    mua = np.array(geom['mua'], float)
    musp = np.array(geom['musp'], float)
    mus = musp / (1.0 - g)               # scattering coefficient for HG stepping
    nlay = len(mua); cortex_z = geom['cortex_z']
    bounds_ext = np.append(bounds, z_max)

    x=np.zeros(N); y=np.zeros(N); z=np.zeros(N)
    ux=np.zeros(N); uy=np.zeros(N); uz=np.ones(N)
    L=np.zeros((N,nlay)); layer=np.zeros(N,int)
    active=np.ones(N,bool); detected=np.zeros(N,bool); exit_r=np.full(N,-1.0)

    for _ in range(max_iter):
        idx=np.where(active)[0]
        if idx.size==0: break
        lay=layer[idx]; zz=z[idx]; uzz=uz[idx]
        step=-np.log(rng.random(idx.size))/mus[lay]
        lower=bounds_ext[lay+1]; upper=bounds[lay]
        target=np.where(uzz>0, lower, upper)
        with np.errstate(divide='ignore',invalid='ignore'):
            dist_b=(target-zz)/uzz
        dist_b=np.where(np.abs(uzz)<1e-12, np.inf, dist_b)
        dist_b=np.where(dist_b<0, np.inf, dist_b)
        hit=step>=dist_b
        s=np.where(hit, dist_b, step)
        L[idx,lay]+=s
        x[idx]+=s*ux[idx]; y[idx]+=s*uy[idx]; z[idx]+=s*uzz
        # boundary hitters
        hb=idx[hit]
        if hb.size:
            lay_hb=layer[hb]; up=uz[hb]<0
            z[hb]=np.where(up, bounds[lay_hb], bounds_ext[lay_hb+1])
            at_surf=up&(lay_hb==0)
            su=hb[at_surf]
            if su.size:
                ci=-uz[su]; si=np.sqrt(np.clip(1-ci**2,0,1)); st_=n_tissue*si
                tir=st_>=1.0; ct=np.sqrt(np.clip(1-np.minimum(st_,1.0)**2,0,1))
                rs=(n_tissue*ci-ct)/(n_tissue*ci+ct+1e-12)
                rp=(n_tissue*ct-ci)/(n_tissue*ct+ci+1e-12)
                R=np.where(tir,1.0,0.5*(rs**2+rp**2))
                refl=rng.random(su.size)<R
                uz[su[refl]]=-uz[su[refl]]
                ex=su[~refl]; detected[ex]=True; active[ex]=False
                exit_r[ex]=np.sqrt(x[ex]**2+y[ex]**2)
            other=hb[~at_surf]
            if other.size:
                lo=layer[other]; upo=uz[other]<0
                layer[other]=np.where(upo, lo-1, lo+1)
                lost=(~upo)&(layer[other]>=nlay)
                if lost.any(): active[other[lost]]=False
                layer[other]=np.clip(layer[other],0,nlay-1)
        # scatterers: HG deflection relative to current direction
        sc=idx[~hit]
        if sc.size:
            xi=rng.random(sc.size)
            if g==0:
                cost=2*xi-1
            else:
                cost=(1+g*g-((1-g*g)/(1-g+2*g*xi))**2)/(2*g)
            cost=np.clip(cost,-1,1); sint=np.sqrt(1-cost*cost)
            phi=2*np.pi*rng.random(sc.size); cph=np.cos(phi); sph=np.sin(phi)
            uxs=ux[sc]; uys=uy[sc]; uzs=uz[sc]
            near=np.abs(uzs)>0.99999
            d=np.sqrt(np.clip(1-uzs*uzs,1e-12,1))
            nux=sint*(uxs*uzs*cph-uys*sph)/d+uxs*cost
            nuy=sint*(uys*uzs*cph+uxs*sph)/d+uys*cost
            nuz=-sint*cph*d+uzs*cost
            # near-vertical fallback
            nux=np.where(near, sint*cph, nux)
            nuy=np.where(near, sint*sph, nuy)
            nuz=np.where(near, np.sign(uzs)*cost, nuz)
            nrm=np.sqrt(nux*nux+nuy*nuy+nuz*nuz)
            ux[sc]=nux/nrm; uy[sc]=nuy/nrm; uz[sc]=nuz/nrm
        # path-length cutoff (weight already negligible)
        Lt=L[idx].sum(axis=1)
        active[idx[Lt>L_max]]=False
    # score
    w=np.exp(-(L*mua[None,:]).sum(axis=1))
    ci=[k for k,b in enumerate(bounds) if b>=cortex_z-1e-9]
    Lcortex=L[:,ci[0]:].sum(axis=1) if ci else np.zeros(N)
    Ltot=L.sum(axis=1)
    det=detected&(exit_r>=0)
    out={}
    for sds in sds_centers:
        m=det&(np.abs(exit_r-sds)<=half_width); ww=w[m]
        out[str(sds)]=dict(sum_wLcortex=float((ww*Lcortex[m]).sum()),
                           sum_wLtot=float((ww*Ltot[m]).sum()),
                           sum_w=float(ww.sum()), count=int(m.sum()))
    out['_meta']=dict(N=N,seed=seed,g=g,detected=int(det.sum()),
                      mean_Ltot=float(Ltot[det].mean()) if det.any() else 0)
    return out

# paper two-layer model at 760 nm (default_adult): 12 mm superficial + cortex
GEOM2L = dict(bounds=[0.0,12.0], mua=[0.0128,0.0192], musp=[1.08,0.82], cortex_z=12.0)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description=(
            "Two-layer Henyey-Greenstein (g=0.9) white Monte Carlo for the cortical "
            "sensitivity fraction f_cortex (and kappa_PV = 1/f_cortex) of the paper's "
            "two-layer head model, evaluated at SDS = 25, 30, 35, 38, 40 mm."),
        epilog=(
            "Examples:\n"
            "  python mc_2layer.py                       # run with the paper defaults\n"
            "  python mc_2layer.py --N 900000 --out f_cortex.json\n"
            "\nThe converged f_cortex is seed-insensitive; the paper used N = 3e5-9e5, "
            "seed = 1, g = 0.9. Results (per-SDS f_cortex, kappa_PV and photon counts) "
            "are written to the --out JSON file."),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-N', '--N', type=lambda v: int(float(v)), default=300000,
                        help="number of photons to launch (default: 300000)")
    parser.add_argument('-s', '--seed', type=int, default=1,
                        help="RNG seed (default: 1; converged result is seed-insensitive)")
    parser.add_argument('-g', '--g', type=float, default=0.9,
                        help="Henyey-Greenstein anisotropy factor (default: 0.9)")
    parser.add_argument('-o', '--out', default='mc_2layer_fcortex.json',
                        help="output JSON file for the per-SDS results "
                             "(default: mc_2layer_fcortex.json)")
    args = parser.parse_args()
    N, seed, g, out = args.N, args.seed, args.g, args.out
    sds=[25,30,35,38,40]; t0=time.time()
    res=run(GEOM2L,N,seed,sds,g=g); res['_meta']['secs']=round(time.time()-t0,1)
    json.dump(res,open(out,'w'))
    print(f"g={g} N={N} seed={seed} detected={res['_meta']['detected']} meanL={res['_meta']['mean_Ltot']:.0f}mm ({res['_meta']['secs']}s)")
    print(f"wrote results to {out}")
    for s in sds:
        d=res[str(s)]; fc=d['sum_wLcortex']/d['sum_wLtot'] if d['sum_wLtot']>0 else float('nan')
        print(f"  SDS={s}: count={d['count']:>6d}  f_cortex={fc:.4f}  kappa_PV={1/fc:.2f}" if fc==fc else f"  SDS={s}: no photons")
