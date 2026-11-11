# SIXEQS_IE_NPHASE verification

These checks call production kernels and inspect actual saved states. They cover
fluid pressure, pressure-temperature and same-substance chemical equilibrium for
EOS 0–3. EOS 2 is the constant-reference affine Mie–Grüneisen model. No separate
thermal-energy correction is built or tested.

Run from the repository root, using a CPU/MPI AMReX checkout:

```sh
make -C exec/ShockBubble-WaterAir-6Eq-Nphase -j4 AMREX_HOME=/path/to/amrex DIFFUSION=true
python3 verification/sixeq_ie_native.py --amrex /path/to/amrex --class-reflux
python3 verification/stiffened_shock_conditions.py --amrex /path/to/amrex
python3 verification/sixeq_ie_applications.py --amrex /path/to/amrex
python3 scripts/test_cases.py --quiet -- AMREX_HOME=/path/to/amrex -j4
```

The shell suite `bash scripts/test_cases.sh` is the other supported entry point.
The native driver reuses only the AMReX library objects of the existing 2D CPU/MPI
build; it recompiles every kernel-test translation unit and any required case
parameter code. `--class-reflux` additionally links the matching two-phase
shock-bubble application objects. A different dimension/backend needs a separate
build and is not covered by these commands.

Use `--run-root /path/to/evidence` on the Python verification drivers to retain
commands, build/run logs and JSON results in a chosen directory. Otherwise they
use temporary directories. The existing suite writes `scripts/logs` and
`scripts/test_summary.csv`. Its `TIMEOUT` status means only that a brief run
survived its time window; it does not establish final-time accuracy or stability.

## What is checked

- Pressure relaxation: 320 EOS/trace cases, energy and volume closure,
  endpoint-work identities, fixed points, tiny/absent phases, positive-residual
  Newton starts, roundoff/domain branches and invalid-input rejection.
- Both HLLC formulations: 492 jump, flux, energy-offset and simplex identities.
- P-T/compatible diffusion: 800 joint-equilibrium projections for each of two and
  three phases, interface traction and heat flux, discrete energy/heating
  consistency and variable viscosity tensors.
- Chemical equilibrium: three/five phases, group-mass and energy conservation,
  chemical-potential/entropy conditions and phase disappearance; parser rejection
  of invalid pairs/reference data and repeat-safe reference-unit conversion.
- Explicit diffusion timestep: periodic heat updates on anisotropic 2D/3D grids,
  preserving the discrete maximum principle below the bound and exposing
  checkerboard growth above it.
- AMR compatible energy source: an actual flux-register/class correction checks
  the traction, work and heat slots against the conservative energy change.
- Shock initialization: 362 Rankine–Hugoniot/branch cases and 180 comparisons with
  the original gamma/pinf interface.
- Applications: 23 short runs of the existing shock-bubble case, including EOS
  0–3, P/P-T, AMR/subcycling, Phase-Field, compatible diffusion and four time
  integrators. Every final checkpoint is audited with the production EOS over
  the composite grid. The auditor currently assumes refinement ratio 2.
- Restart: two AMR runs compared at the same physical time, with and without
  subcycling; saved-state checks and normalized state difference below `1e-8`.
  Adaptive timestep estimates are reconstructed at restart, so bitwise identity
  is not promised.

## Advection admissibility regression

```sh
python3 verification/sixeq_ie_advection.py --amrex /path/to/amrex
# Optional two-process check of a selected variant:
python3 verification/sixeq_ie_advection.py --amrex /path/to/amrex --ranks 2 --only WENO5B
```

This driver uses the existing five-phase advection case, runs Godunov, MUSCL,
WENO3B, WENO5 and WENO5B to the same physical time (`0.01` by default), and audits
initial/final composite-grid checkpoints. It checks every partial mass, total
energy, phase volumes, bounds, common pressure and uniform velocity. Native
reconstruction tests exercise EOS 0–3, tiny/exactly absent phases, the interior
quadrature state, stretched pure liquid and smooth-grid refinement.

The five previously failing variants are addressed by three explicit changes:

- This benchmark uses `run.do_subcycle=0`: fine/coarse fluxes agree at every RK
  stage. Both refinement levels remain. The general subcycled reflux algorithm
  is unchanged; selecting subcycling again can reproduce its positivity failure.
- WENO3B/WENO5B reconstruct intrinsic phasic pressure and form `alpha*p` with the
  limited fraction. WENO5 also uses a shared convex reconstruction limiter.
  Fractions/masses and the EOS thermal domain are bounded at both faces and the
  interior quadrature point. Limiting acts on reconstructed polynomials; evolved
  cell averages are not clipped, floored or renormalized after an update.
- IE WENO3B/WENO5/WENO5B obey the unsplit quadrature timestep cap
  `dt <= 0.9*w/(cmax*sum_d(1/dx_d))`, with `w=1/6` for WENO3B and `w=1/12` for
  fifth-order WENO. The user CFL remains an upper bound. This cap alone is not a
  positivity proof for every nonlinear multiphysics update; the state checks stay on.

The component-wise WENO5 and shared-weight WENO5B reconstructions remain distinct.
Both may reduce their reconstruction near a phase boundary. Neither adds the
separate thermal-energy correction algorithm. An admissible extension of the
unobservable intrinsic pressure of an absent phase leaves its weighted pressure,
mass and energy zero, including a stretched liquid with an absent gas.

The upstream subcycling issue was reproduced at commit `408698a`: even Godunov
reached a transient `alpha=-5.038754462e-6`, although its saved end-of-step states
were positive. Final-state checks alone would miss that issue. The regression
therefore keeps the stage admissibility checks enabled.

The existing CI test workflow is disabled (`.github/workflows/compas-test.yml.disabled`)
and pins AMReX 26.01; local validation used AMReX 26.02-53 and CPU/MPI. GPU backends,
the pinned CI release and a full-resolution complete advection period have not
been validated by these short regressions.
