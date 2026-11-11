#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

from pathlib import Path
import subprocess,json,os,re,argparse,tempfile
ap=argparse.ArgumentParser(description="Short existing shock-bubble IE runs with production-EOS checkpoint audits")
ap.add_argument('--amrex',type=Path,required=True)
ap.add_argument('--run-root',type=Path)
args=ap.parse_args()
R=Path(__file__).resolve().parents[1]
V=(args.run_root or Path(tempfile.mkdtemp(prefix='compas-ie-apps-'))).resolve();V.mkdir(parents=True,exist_ok=True)
case=R/'exec/ShockBubble-WaterAir-6Eq-Nphase';cache=case/'tmp_build_dir';amrex=args.amrex.resolve();
with (V/'application-build.log').open('w') as log:
 subprocess.run(['make','-j4','AMREX_HOME='+str(amrex),'DIFFUSION=true'],cwd=case,stdout=log,stderr=subprocess.STDOUT,check=True)
env=dict(os.environ,OMPI_MCA_btl='self')
includes=[cache/'s/2d.gnu.MPI.EXE',case/'prob',R/'include']+[p for p in (R/'include').iterdir() if p.is_dir()]+[amrex/'Src'/p for p in ('Base','Base/Parser','Boundary','AmrCore','LinearSolvers/MLMG')]
objects=[p for p in (cache/'o/2d.gnu.MPI.EXE').glob('*.o') if p.name.startswith(('AMReX','amrex_'))]
cmd=['mpicxx','-std=c++17','-O2','-DOMPI_SKIP_MPICXX','-DPHYSICS=SIXEQS_IE_NPHASE','-DNPHASE=2','-DDIFFUSION=true','-DADVECTION=true','-DNONCONSERVATIVE=true']+['-I'+str(p) for p in includes]+[str(R/'verification/sixeq_ie_testcase_audit.cpp'),str(R/'source/Parm.cpp')]+[str(p) for p in objects]+['-o',str(V/'audit')]
with (V/'audit-build.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
results={}
configs=[(eos,mode,0,0,0,'TVD-RK3',0) for eos in range(4) for mode in ['Pressure','PressureTemperature']]
configs += [(1,mode,1,sub,pf,ti,compatible) for mode,compatible in [('Pressure',0),('PressureTemperature',0),('PressureTemperature',1)] for sub in (0,1) for pf in (0,1) for ti in ['TVD-RK3']]
configs += [(1,'PressureTemperature',1,0,0,ti,1) for ti in ['ForwardEuler','TVD-RK2','RK4']]
for eos,mode,level,sub,pf,ti,compatible in configs:
 name=f'eos{eos}-{mode}-L{level}-sub{sub}-pf{pf}-{ti}-compatible{compatible}';d=V/name;d.mkdir(exist_ok=True)
 opts=[f'prob.EOS.EOS={eos}',f'Physics.Relaxation={mode}',f'Physics.CompatibleDiffusion={compatible}','max_step=10','stop_time=1',f'amr.max_level={level}','amr.n_cell=32 16 8','amr.max_grid_size=16','amr.blocking_factor_x=8','amr.blocking_factor_y=8',f'run.do_subcycle={sub}',f'PhaseField.phase_field={pf}',f'run.TimeIntegrator={ti}','run.cfl=.1','amr.regrid_int=2','amr.plot_int=-1','amr.t_write_interval=0','amr.chk_int=5','amr.chk_file='+str(d/'chk'),'amr.case_name='+str(d/'plot/case'),'run.user_output_int=0','amrex.signal_handling=0']
 if eos==2:opts+=['prob.EOS.gGamma=3.4 .4','prob.EOS.pref=-2.64e9 0','prob.EOS.eref=1000 -2000']
 if eos==3:opts+=['prob.EOS.b=.0001 .00001','prob.EOS.q=1000 -2000']
 if compatible:opts+=['prob.EOS.mu=.001 .000018','prob.EOS.kappa=.6 .026']
 cmd=[str(case/'main2d.gnu.MPI.ex'),str(case/'prob/inputs')]+opts
 (d/'command.json').write_text(json.dumps(cmd,indent=2))
 p=subprocess.run(cmd,cwd=d,env=env,capture_output=True,text=True,timeout=120);out=p.stdout+p.stderr;(d/'run.log').write_text(out)
 result={'run_ok':p.returncode==0}
 if p.returncode:
  result['diagnostic']=next((l for l in out.splitlines() if 'amrex::Abort' in l),out[-800:])
 else:
  chk=sorted(p for p in d.glob('chk*') if re.fullmatch(r'chk[0-9]+',p.name))[-1];finest=int((chk/'Header').read_text().splitlines()[1])
  cmd=[str(V/'audit'),str(case/'prob/inputs')]+opts+[f'audit.checkpoint={chk}',f'audit.levels={finest+1}']
  p=subprocess.run(cmd,cwd=d,env=env,capture_output=True,text=True,timeout=60);out=p.stdout+p.stderr;(d/'audit.log').write_text(out);result['audit_ok']=p.returncode==0
  m=re.search(r'CASE_AUDIT (\{.*\})',out)
  if m:result['audit']=json.loads(m[1])
  if p.returncode:result['diagnostic']=out[-1000:]
 results[name]=result;(V/'results.json').write_text(json.dumps(results,indent=2)+'\n');print('PASS' if result.get('audit_ok') else 'FAIL',name,flush=True)

# Restart from the intermediate checkpoint of a matched continuous run.
for sub in (0,1):
 name=f'eos1-PressureTemperature-L1-sub{sub}-pf0-TVD-RK3-compatible1'
 original=V/name;restart=V/('restart-sub'+str(sub));restart.mkdir(exist_ok=True)
 cmd=json.loads((original/'command.json').read_text())
 target_time=float((original/'chk00010/Header').read_text().splitlines()[4].split()[0])
 cmd=[v.replace(str(original),str(restart)) for v in cmd if not v.startswith(('stop_time=','max_step='))]+['amr.restart='+str(original/'chk00005'),'stop_time='+repr(target_time),'max_step=30','amr.chk_int=1']
 p=subprocess.run(cmd,cwd=restart,env=env,capture_output=True,text=True,timeout=120)
 (restart/'run.log').write_text(p.stdout+p.stderr)
 result={'run_ok':p.returncode==0}
 if not p.returncode:
  comparison=V/'checkpoint-compare'
  command=['mpicxx','-std=c++17','-O2','-DOMPI_SKIP_MPICXX']+['-I'+str(p) for p in includes]+[str(R/'verification/sixeq_ie_checkpoint_compare.cpp')]+[str(p) for p in objects]+['-o',str(comparison)]
  with (V/'compare-build.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
  restart_final=sorted(p for p in restart.glob('chk*') if re.fullmatch(r'chk[0-9]+',p.name))[-1]
  original_header=(original/'chk00010/Header').read_text().splitlines();restart_header=(restart_final/'Header').read_text().splitlines()
  result['max_time_error']=max(abs(float(a)-float(b)) for a,b in zip(original_header[4].split(),restart_header[4].split()))
  p=subprocess.run([str(comparison),'left='+str(original/'chk00010'),'right='+str(restart_final),'levels=2','amrex.signal_handling=0'],cwd=restart,env=env,capture_output=True,text=True,timeout=60)
  out=p.stdout+p.stderr;(restart/'comparison.log').write_text(out)
  m=re.search(r'CHECKPOINT_COMPARE (\{.*\})',out)
  result['comparison']=json.loads(m[1]) if m else {}
  result['audit_ok']=p.returncode==0 and result['comparison'].get('max_relative_error',1)<1e-8 and result['max_time_error']<1e-12
  audit_cmd=[str(V/'audit')]+cmd[1:]+['audit.checkpoint='+str(restart_final),'audit.levels=2']
  p=subprocess.run(audit_cmd,cwd=restart,env=env,capture_output=True,text=True,timeout=60)
  (restart/'audit.log').write_text(p.stdout+p.stderr)
  result['state_audit_ok']=p.returncode==0
  result['audit_ok']=result['audit_ok'] and result['state_audit_ok']
 results['restart-sub'+str(sub)]=result;(V/'results.json').write_text(json.dumps(results,indent=2)+'\n');print('PASS' if result.get('audit_ok') else 'FAIL','restart-sub'+str(sub),flush=True)
if any(not r.get('audit_ok') for r in results.values()):raise SystemExit('One or more application/audit checks failed; see results.json')
print('Evidence:',V)
