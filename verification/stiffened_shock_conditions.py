# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Check the shared shock initializer against RH jumps and its legacy interface.

Build an existing CPU testcase first. This compiles a separate verification
executable using that testcase's AMReX library objects; no inputs are edited.
"""
from pathlib import Path
import argparse,json,os,subprocess,tempfile

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,default=Path('exec/ShockBubble-WaterAir-6Eq-Nphase'))
    parser.add_argument('--run-root',type=Path)
    parser.add_argument('--amrex',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    case=args.case if args.case.is_absolute() else root/args.case
    cache=case/'tmp_build_dir';suffix='2d.gnu.MPI.EXE'
    amrex=args.amrex or root.parent/'amrex'
    run=args.run_root or Path(tempfile.mkdtemp(prefix='compas-shock-check-'))
    run=run.resolve();run.mkdir(parents=True,exist_ok=True)
    objects=[p for p in (cache/'o'/suffix).glob('*.o') if p.name.startswith(('AMReX','amrex_'))]
    if not objects:raise RuntimeError('Build the selected 2D CPU testcase first')
    includes=[cache/'s'/suffix,root/'include']+[amrex/'Src'/p for p in ['Base','Base/Parser','Boundary','AmrCore','LinearSolvers/MLMG']]
    command=['mpicxx','-std=c++17','-O2','-DOMPI_SKIP_MPICXX']+['-I'+str(p) for p in includes]+[str(root/'verification/stiffened_shock_conditions.cpp')]+[str(p) for p in objects]+['-o',str(run/'shock-check')]
    (run/'build-command.json').write_text(json.dumps(command,indent=2)+'\n')
    with (run/'build.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    p=subprocess.run([str(run/'shock-check'),'amrex.signal_handling=0'],cwd=run,env=dict(os.environ,OMPI_MCA_btl='self'),capture_output=True,text=True)
    (run/'run.log').write_text(p.stdout+p.stderr)
    if p.returncode:raise RuntimeError(p.stdout+p.stderr)
    print(p.stdout);print('Evidence:',run)

if __name__=='__main__':main()
