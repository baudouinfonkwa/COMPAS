#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
"""Compile fluid IE kernel tests against a built 2D CPU/MPI AMReX cache.

Build exec/ShockBubble-WaterAir-6Eq-Nphase first with DIFFUSION=true.
Only AMReX library objects are reused; every application-dependent object is rebuilt.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--amrex', type=Path, required=True)
    ap.add_argument('--case', type=Path, default=Path('exec/ShockBubble-WaterAir-6Eq-Nphase'))
    ap.add_argument('--suffix', default='2d.gnu.MPI.EXE')
    ap.add_argument('--class-reflux', action='store_true', help='Also link the built two-phase shock-bubble classes to test AMR source slots')
    ap.add_argument('--run-root', type=Path)
    args = ap.parse_args()
    args.amrex = args.amrex.resolve()
    case = (ROOT / args.case).resolve()
    cache = case / 'tmp_build_dir'
    run = (args.run_root or Path(tempfile.mkdtemp(prefix='compas-ie-native-'))).resolve()
    run.mkdir(parents=True, exist_ok=True)
    objects = sorted(p for p in (cache / 'o' / args.suffix).glob('*.o')
                     if p.name.startswith(('AMReX', 'amrex_')))
    if not objects:
        raise RuntimeError('Build the selected 2D CPU/MPI case first')
    env = dict(os.environ, OMPI_MCA_btl='self')
    results = {}

    def compile_test(name, nphase, diffusion='true', testcase=case, parm=False, application_objects=False):
        includes = [cache / 's' / args.suffix, testcase / 'prob', ROOT / 'include']
        includes += [p for p in (ROOT / 'include').iterdir() if p.is_dir()]
        includes += [args.amrex / 'Src' / p for p in
                     ('Base', 'Base/Parser', 'Boundary', 'AmrCore', 'LinearSolvers/MLMG')]
        exe = run / f'{name}-{nphase}-{diffusion}'
        command = ['mpicxx', '-std=c++17', '-O2', '-DOMPI_SKIP_MPICXX',
                   '-DPHYSICS=SIXEQS_IE_NPHASE', f'-DNPHASE={nphase}',
                   '-DADVECTION=true', '-DNONCONSERVATIVE=true', '-DDIFFUSION=' + diffusion]
        command += ['-I' + str(p) for p in includes]
        command += [str(ROOT / 'verification' / (name + '.cpp'))]
        if parm:
            command += [str(ROOT / 'source/Parm.cpp')]
        linked = [p for p in (cache / 'o' / args.suffix).glob('*.o') if p.name != 'main.o'] if application_objects else objects
        command += [str(p) for p in linked] + ['-o', str(exe)]
        (run / (exe.name + '.command.json')).write_text(json.dumps(command, indent=2))
        with (run / (exe.name + '.build.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        return exe

    def execute(exe, label, options=(), expected=None, marker=None):
        p = subprocess.run([str(exe), *options, 'amrex.signal_handling=0'], cwd=run,
                           env=env, capture_output=True, text=True, timeout=60)
        output = p.stdout + p.stderr
        (run / (label + '.log')).write_text(output)
        if expected is not None:
            if p.returncode == 0 or expected not in output:
                raise RuntimeError(f'{label}: expected rejection {expected!r}\n{output[-4000:]}')
            results[label] = {'rejected': True, 'diagnostic': expected}
        else:
            if p.returncode:
                raise RuntimeError(f'{label}: failed\n{output[-4000:]}')
            match = re.search(re.escape(marker or '') + r' (\{.*\})', output) if marker else None
            results[label] = json.loads(match[1]) if match else {'passed': True}
        (run / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        print('PASS', label, flush=True)

    exe = compile_test('sixeq_ie_general_relaxation', 3)
    execute(exe, 'pressure-eos0-3', marker='RELAXATION_RESULTS')
    for mode, diagnostic in [('energy', 'produced an inadmissible state'),
                             ('temperature', 'incoming phase is below its EOS temperature bound'),
                             ('fraction', 'invalid incoming volume fraction or partial mass'),
                             ('covolume', 'covolume bound violated')]:
        execute(exe, 'reject-' + mode, ['invalid=' + mode], expected=diagnostic)
    exe = compile_test('sixeq_ie_hllc_variants', 3)
    execute(exe, 'hllc-eos0-3', marker='HLLC_RESULTS')
    for choice in ('ConservativePart', 'Consistent', '0', '1'):
        execute(exe, 'hllc-parse-' + choice, ['audit_mode=parse', 'Physics.HLLCImplementation=' + choice])
    execute(exe, 'hllc-reject', ['audit_mode=parse', 'Physics.HLLCImplementation=invalid'],
            expected='Physics.HLLCImplementation must be')
    for nphase in (3, 5):
        exe = compile_test('sixeq_ie_chemical', nphase, testcase=ROOT / 'exec/_Tests/Diffusion')
        execute(exe, 'chemical-' + str(nphase), marker='CHEMICAL_RESULTS')
    for nphase, deck in [(2, 'inputs.Couette2Layer-6Eq-N'), (3, 'inputs.Advection-Viscous-6Eq-N')]:
        testcase = ROOT / 'exec/_Tests/Diffusion'
        exe = compile_test('sixeq_ie_compatible_diffusion', nphase, testcase=testcase, parm=True)
        execute(exe, 'diffusion-pt-' + str(nphase), [str(testcase / 'prob' / deck)], marker='COMPATIBLE_DIFFUSION')
    exe = compile_test('sixeq_ie_chemical_flags', 3)
    chemical = ['Physics.pressure_relaxation=1', 'Physics.Relaxation=PressureTemperature',
                'Physics.ChemicalRelaxation=1', 'Physics.ChemicalRelaxation.Pairs=0 1',
                'prob.EOS.Substance=water water air', 'prob.EOS.sref=100 200 300',
                'prob.EOS.T_ref=300', 'prob.EOS.p_ref=100000',
                'test.pressure_unit=100000', 'test.energy_unit=5000']
    execute(exe, 'chemical-reference-units', chemical, marker='CHEMICAL_FLAGS')
    if results['chemical-reference-units']['p_ref'] != 1 or abs(results['chemical-reference-units']['sref0'] - .02) > 1e-15:
        raise RuntimeError('Chemical reference units were not converted correctly')
    invalid = [('pairs-missing', 'Physics.ChemicalRelaxation.Pairs', None, 'Pairs requires'),
               ('pairs-odd', 'Physics.ChemicalRelaxation.Pairs', '0', 'Pairs requires'),
               ('pairs-self', 'Physics.ChemicalRelaxation.Pairs', '0 0', 'invalid phase index or self-pair'),
               ('pairs-range', 'Physics.ChemicalRelaxation.Pairs', '0 3', 'invalid phase index or self-pair'),
               ('pairs-duplicate', 'Physics.ChemicalRelaxation.Pairs', '0 1 1 0', 'duplicate pair'),
               ('substance', 'prob.EOS.Substance', 'water air air', 'same Substance'),
               ('entropy-missing', 'prob.EOS.sref', None, 'requires explicit prob.EOS.sref'),
               ('entropy-short', 'prob.EOS.sref', '100 200', 'requires explicit prob.EOS.sref'),
               ('Tref', 'prob.EOS.T_ref', '0', 'finite positive EOS T_ref'),
               ('pref', 'prob.EOS.p_ref', '0', 'finite positive EOS T_ref'),
               ('tolerance', 'Physics.ChemicalRelaxation.Tolerance', '0', '0<Tolerance'),
               ('sweeps', 'Physics.ChemicalRelaxation.MaxSweeps', '0', 'MaxSweeps>=1'),
               ('pressure-only', 'Physics.Relaxation', 'Pressure', 'requires Physics.Relaxation=PressureTemperature')]
    for label, key, value, diagnostic in invalid:
        options = [v for v in chemical if not v.startswith(key + '=')]
        if value is not None:
            options.append(key + '=' + value)
        execute(exe, 'chemical-reject-' + label, options, expected=diagnostic)
    for diffusion in ('true', 'false'):
        exe = compile_test('sixeq_ie_diffusion_flag', 2, diffusion=diffusion)
        for value, pressure in [(0, 1), (1, 1), (-1, 1), (2, 1), (1, 0)]:
            valid = value == 0 or (value == 1 and pressure == 1 and diffusion == 'true')
            execute(exe, f'diffusion-{diffusion}-flag{value}-p{pressure}',
                    [f'Physics.CompatibleDiffusion={value}', f'Physics.pressure_relaxation={pressure}'],
                    expected=None if valid else 'Physics.CompatibleDiffusion')
    exe = compile_test('sixeq_ie_diffusion_timestep', 2)
    execute(exe, 'diffusion-timestep', marker='DIFFUSION_TIMESTEP')
    for mode, diagnostic in [('coefficient', 'invalid diffusion timestep coefficient'),
                             ('safety', 'requires 0 < run.vnn <= 1')]:
        execute(exe, 'diffusion-timestep-reject-' + mode, ['invalid=' + mode], expected=diagnostic)
    if args.class_reflux:
        if case != ROOT / 'exec/ShockBubble-WaterAir-6Eq-Nphase':
            raise RuntimeError('--class-reflux requires the default two-phase shock-bubble build')
        config = (cache / 's' / args.suffix / 'AMReX_Config.H').read_text()
        if '#define NPHASE 2' not in config or '#define DIFFUSION true' not in config:
            raise RuntimeError('Build the default shock-bubble case with DIFFUSION=true first')
        exe = compile_test('sixeq_ie_compatible_reflux', 2, application_objects=True)
        execute(exe, 'class-reflux', [str(case / 'prob/inputs'), 'amr.n_cell=8 8',
                'amr.max_level=1', 'amr.blocking_factor_x=4', 'amr.blocking_factor_y=4',
                'Physics.CompatibleDiffusion=1', 'Physics.Relaxation=PressureTemperature',
                'amr.case_name=' + str(run / 'plot/case')], marker='COMPATIBLE_REFLUX')
    print('Evidence:', run)


if __name__ == '__main__':
    main()
