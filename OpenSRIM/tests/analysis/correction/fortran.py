from pathlib import Path
import subprocess, shutil
import numpy as np
root=Path(__file__).resolve().parents[4]; im=root/'imsil-omp'; out=Path('/tmp/correction-study'); out.mkdir(exist_ok=True)
sources=[]
for line in (im/'code/srcs.mk').read_text().splitlines():
 s=line.strip().rstrip('\\').strip()
 if s.startswith('${SRCDIR}/') and not s.endswith('/main.f90'):sources.append(im/s.replace('${SRCDIR}','code'))
driver=(root/'OpenSRIM/tests/fortran/generate_extended_scattering_references.f90').read_text()
driver=driver[:driver.index('   ! Keep these input grids')]+'''   open(newunit=unit,file='results.txt',status='replace')
   do ie=1,7
      e=10._REAL8**(ie-5)
      do ip=1,215
         if(ip<=201) then
            p=0.01_REAL8+(ip-1)*(30._REAL8-0.01_REAL8)/200
         else
            p=rmax*(1-10._REAL8**(-(ip-200)))
         end if
         call cmAngle(e,p,1,2,theta,method='LegendreChiCorr',nabsc=4)
         write(unit,'(3es25.16)') e,p,theta
      end do
   end do
   close(unit)
end program generate_extended_scattering_references
'''
raw=(im/'code/sn/cmscatter.f90').read_text(); expr='FOUR*interaction/(base*h*(h+p_save*base))'
for mode,replacement in [('current',expr),('zero','ZERO'),('subtract','FOUR*(ONE/base-p_save/h)')]:
 stage=out/mode;stage.mkdir(exist_ok=True);(stage/'mod').mkdir(exist_ok=True)
 for name in ('gauss_legendre.dat','rhoe.dat','rhozbl.dat','scattab','screentab'):
  if not (stage/name).exists():(stage/name).symlink_to(im/'tables'/name)
 (stage/'screencoefs').mkdir(exist_ok=True)
 (stage/'screencoefs/NLHLIN_python_coeffs.dat').write_text('# nexp=3\n'+(root/'OpenSRIM/data/NLHlin/NLHlin_3eV.dat').read_text())
 (stage/'cmscatter.f90').write_text(raw.replace(expr,replacement));(stage/'driver.f90').write_text(driver)
 src=[stage/'cmscatter.f90' if s.name=='cmscatter.f90' else s for s in sources]
 with (stage/'build.log').open('w') as log:
  subprocess.run(['gfortran','-O0','-J',str(stage/'mod'),'-cpp','-o',str(stage/'run'),*map(str,src),str(stage/'driver.f90')],stdout=log,stderr=log,check=True)
 for z1,z2 in [(5,14),(1,1),(79,79),(1,79)]:
  with (stage/'run.log').open('w') as log:subprocess.run([str(stage/'run'),'NLHlin',str(z1),str(z2),'512'],cwd=stage,stdout=log,stderr=log,check=True)
  shutil.copy(stage/'results.txt',stage/f'{z1}_{z2}.txt')
 print(mode,'done',flush=True)
for mode in ['zero','subtract']:
 for group in ['regular','boundary']:
  vals=[]
  for pair in ['5_14','1_1','79_79','1_79']:
   a=np.loadtxt(out/'current'/f'{pair}.txt');b=np.loadtxt(out/mode/f'{pair}.txt')
   for i,(x,y) in enumerate(zip(a,b)):
    if (i%215<201)==(group=='regular') and np.isfinite(x[2]) and np.isfinite(y[2]) and x[2]!=0:vals.append((abs(x[2]-y[2]),abs((x[2]-y[2])/x[2]),pair,*x,y[2]))
  print(mode,group,'abs',max(vals,key=lambda x:x[0]),'rel',max(vals,key=lambda x:x[1]))
