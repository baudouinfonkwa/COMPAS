#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
"""Check five-phase IE reconstruction and equal-time AMR advection conservation."""
from pathlib import Path
import argparse,csv,json,os,re,subprocess,tempfile


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--amrex',type=Path,required=True)
    ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--exe',type=Path)
    ap.add_argument('--cache',type=Path)
    ap.add_argument('--native-source',type=Path)
    ap.add_argument('--run-root',type=Path)
    ap.add_argument('--stop-time',type=float,default=.01)
    ap.add_argument('--ranks',type=int,default=1)
    ap.add_argument('--only',choices=['GODUNOV','MUSCL','WENO3B','WENO5','WENO5B'])
    args=ap.parse_args();root=args.repo.resolve();case=root/'exec/_Tests/Advection';amrex=args.amrex.resolve()
    run=(args.run_root or Path(tempfile.mkdtemp(prefix='compas-ie-advection-'))).resolve();run.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OMPI_MCA_btl='self' if args.ranks==1 else 'self,sm')
    if args.exe is None:
        with (run/'application-build.log').open('w') as log:
            subprocess.run(['make','-j4','BUILD=6Eq-N','AMREX_HOME='+str(amrex)],cwd=case,stdout=log,stderr=subprocess.STDOUT,check=True)
    exe=(args.exe or case/'main2d.gnu.MPI.6Eq-N.ex').resolve()
    cache=(args.cache or case/'tmp_build_dir').resolve()
    suffix='2d.gnu.MPI.6Eq-N.EXE' if (cache/'s/2d.gnu.MPI.6Eq-N.EXE').is_dir() else '2d.gnu.MPI.EXE'
    config=(cache/'s'/suffix/'AMReX_Config.H').read_text()
    if '#define AMREX_SPACEDIM 2' not in config or '#define NPHASE 5' not in config:
        raise RuntimeError('The cache must be from the five-phase 2D IE advection build')
    includes=[cache/'s'/suffix,case/'prob',root/'include']+[p for p in (root/'include').iterdir() if p.is_dir()]
    includes += [amrex/'Src'/p for p in ['Base','Base/Parser','Boundary','AmrCore','LinearSolvers/MLMG']]
    objects=[p for p in (cache/'o'/suffix).glob('*.o') if p.name.startswith(('AMReX','amrex_'))]
    if not objects:raise RuntimeError('Build the selected 2D CPU/MPI five-phase IE case first')
    base=['mpicxx','-std=c++17','-O2','-DOMPI_SKIP_MPICXX','-DPHYSICS=SIXEQS_IE_NPHASE','-DNPHASE=5','-DDIFFUSION=false','-DADVECTION=true','-DNONCONSERVATIVE=true']+['-I'+str(p) for p in includes]
    def compile_test(source,name,parm=False):
        cmd=base+[str(source)]+([str(root/'source/Parm.cpp')] if parm else [])+[str(p) for p in objects]+['-o',str(run/name)]
        with (run/(name+'-build.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        return run/name
    native=compile_test(args.native_source or root/'verification/sixeq_ie_reconstruction.cpp','reconstruction')
    audit=compile_test(root/'verification/sixeq_ie_testcase_audit.cpp','audit',True)
    def execute(cmd,d,label,marker=None):
        (d/(label+'-command.json')).write_text(json.dumps([str(v) for v in cmd],indent=2)+'\n')
        p=subprocess.run(cmd,cwd=d,env=env,capture_output=True,text=True,timeout=300)
        out=p.stdout+p.stderr;(d/(label+'.log')).write_text(out)
        if p.returncode:raise RuntimeError(label+' failed: '+out[-2000:])
        if marker:
            m=re.search(re.escape(marker)+r' (\{.*\})',out)
            if not m:raise RuntimeError('Missing '+marker)
            return json.loads(m[1])
    results={'reconstruction':execute([str(native),'amrex.signal_handling=0'],run,'native','RECONSTRUCTION'),'cases':{}}
    for scheme in ['GODUNOV','MUSCL','WENO3B','WENO5','WENO5B']:
        if args.only and scheme!=args.only:continue
        d=run/scheme;d.mkdir(exist_ok=True)
        opts=['run.do_subcycle=1','max_step=10000','stop_time='+str(args.stop_time),'FiniteVolume.Scheme='+scheme,
              'amr.plot_int=1000000','amr.t_write_interval=0','run.output_vars=a1 P','amr.chk_int=1000000',
              'run.user_output_int=1','amr.chk_file='+str(d/'chk'),'amr.case_name='+str(d/'plot/case'),
              'amrex.signal_handling=0']
        deck=str(case/'prob/inputs.Advection-6Eq-N')
        execute((['mpirun','-n',str(args.ranks)] if args.ranks>1 else [])+[str(exe),deck]+opts,d,'run')
        checkpoints=sorted(p for p in d.glob('chk*') if re.fullmatch(r'chk[0-9]+',p.name))
        if len(checkpoints)<2:raise RuntimeError('Missing initial/final checkpoints for '+scheme)
        data=[]
        for n,checkpoint in enumerate([checkpoints[0],checkpoints[-1]]):
            header=(checkpoint/'Header').read_text().splitlines();levels=int(header[1])+1
            options=opts+['audit.checkpoint='+str(checkpoint),'audit.levels='+str(levels),
                          'audit.expected_pressure=1','audit.expected_velocity_x=1']
            if not n:options+=['audit.check_initial_floors=1']
            data.append(execute([str(audit),deck]+options,d,'audit-'+str(n),'CASE_AUDIT'))
        initial,final=data
        time=float((checkpoints[-1]/'Header').read_text().splitlines()[4].split()[0])
        mass_error=max(abs(a-b)/max(1.,abs(a)) for a,b in zip(initial['phase_mass_sums'],final['phase_mass_sums']))
        energy_error=abs(final['total_energy_sum']-initial['total_energy_sum'])/max(1.,abs(initial['total_energy_sum']))
        rows=list(csv.DictReader((d/'plot/case/advection.csv').open()));row=rows[-1]
        volume_error=max(abs(float(row[k])-float(rows[0][k])) for k in row if k.startswith('V_'))
        result={'time':time,'steps':int(row['step']),'max_mass_drift':mass_error,'total_energy_drift':energy_error,
                'max_phase_volume_drift':volume_error,'pressure_error':float(row['dP']),'velocity_error':float(row['du']),
                'final_state':final}
        if abs(time-args.stop_time)>1e-12 or mass_error>1e-11 or energy_error>1e-11 or volume_error>1e-11 or final['minimum_alpha']<0 or final['minimum_mass']<0 or result['pressure_error']>1e-10 or result['velocity_error']>1e-10:
            raise RuntimeError(scheme+' failed equal-time conservation/equilibrium checks: '+str(result))
        results['cases'][scheme]=result;(run/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        print('PASS',scheme,'steps',result['steps'],'mass',mass_error,'energy',energy_error,flush=True)
    print('Evidence:',run)


if __name__=='__main__':main()
