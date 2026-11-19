import sys, types, json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
Path('/tmp/correction-study').mkdir(exist_ok=True)
sys.path[:0]=[str(ROOT),str(ROOT/'OpenSRIM/tests')]
from reference_data import python_modules
from scattering_support import numba_coefficients
from OpenSRIM.simulators.opentrim import cm_scatter as nb
mods=python_modules(); py=mods['cm_scatter']

def clone(module, source, suffix):
    m=types.ModuleType(module.__name__+'_'+suffix); m.__package__=module.__package__
    sys.modules[m.__name__]=m
    exec(compile(source, module.__file__, 'exec'),m.__dict__)
    return m
ps=Path(py.__file__).read_text(); ns=Path(nb.__file__).read_text()
variants={}
for name,mod,src in [('Python',py,ps),('Numba',nb,ns)]:
    expr='rho*chi / (base*g*(g+base))' if name=='Python' else '(r0/e*chi) / (base*h*(h+p*base))'
    direct='1/base - 1/g' if name=='Python' else '1/base - p/h'
    variants[name]=[mod,clone(mod,src.replace('correction = '+expr,'correction = np.zeros_like(u)'), 'zero'),clone(mod,src.replace('correction = '+expr,'correction = '+direct),'subtract')]
rows=[]
for pair in [(5,14),(1,1),(79,79),(1,79)]:
    norm,coefs=numba_coefficients('NLHlin',*pair)
    screen=mods['nlhlin'].NLHlin_screen(*pair,norm); R=screen.rmax
    # Regular test range plus logarithmic and cutoff stress probes.
    grids={'regular':np.linspace(.01,30,201),'boundary':R*(1-10.**(-np.arange(2,16,dtype=float)))}
    for n in [4,32]:
        for m in variants['Python']:m.setup(n)
        for group,points in grids.items():
            for e in np.logspace(-4,2,7):
                for p in points:
                    for backend,ms in variants.items():
                        vals=[m.scatter_integrals(e,p,screen) if backend=='Python' else m.scatter_integrals(e,p,'NLHlin',coefs,n) for m in ms]
                        theta=[v[0] if backend=='Python' else np.pi-v[0] for v in vals]
                        rows.append([*pair,n,group,float(e),float(p),backend,*theta,max(abs(v[1]-vals[0][1]) for v in vals)])
                    # Correct pre-rationalization expression and historical wrong prefactor.
                    if p<R:
                        r0,_=screen.apsis(e,p); umax=np.sqrt((R-r0)/R)
                        if umax>0:
                            u,w=py.ROOTS_LEGENDRE;u=.5*umax*(u+1);w=.5*umax*w
                            _,chi=py.calc_phi_chi(u,r0,screen)
                            Q=np.sum(w/np.sqrt(2-u*u+r0/(e*p*p)*chi))/np.sum(w/np.sqrt(2-u*u))
                            A=2*np.arccos(p/R);B=2*np.arccos(r0/R)
                            ref=variants['Python'][0].scatter_integrals(e,p,screen)[0]
                            rows.append([*pair,n,group,float(e),float(p),'Formula',ref,A-B*Q,A*(1-Q),0.])
Path('/tmp/correction-study/results.json').write_text(json.dumps(rows))
for backend in ['Python','Numba','Formula']:
 for group in ['regular','boundary']:
  for n in [4,32]:
   arr=[r for r in rows if r[6]==backend and r[3]==group and r[2]==n]
   for col in [8,9]:
    valid=[r for r in arr if np.isfinite(r[7]) and np.isfinite(r[col]) and r[7]!=0]
    ab=max(valid,key=lambda r:abs(r[col]-r[7])); rel=max(valid,key=lambda r:abs((r[col]-r[7])/r[7]))
    print(backend,group,n,col,'MAXABS',abs(ab[col]-ab[7]),ab,'MAXREL',abs((rel[col]-rel[7])/rel[7]),rel,flush=True)
print('max tau difference',max(r[10] for r in rows))

print('nonfinite rows',sum(not np.isfinite(r[7:]).all() for r in rows))
