# Two off-hooked trailers -- MATLAB study

Stabilisation of a TurtleBot3 with **two off-hooked trailers** by the sampled
oscillating controller of Zuyev & Grushkovskaya (arXiv:2403.03341), designed on the
nilpotent approximation, tested against the exact kinematics and the robot's
velocity limits.

## Quick start

```matlab
selftest            % 30 s -- does everything still compute what it did before?
plots               % one simulation: states, controls vs. limits, convergence
compare_physical    % the SAME controller on the exact kinematics (the lab preview)
```

All numbers live in **`params.m`**. Change them there and nowhere else.

## What to run for what

| script | question it answers |
|---|---|
| `selftest` | Did my edit break anything? (run after every change) |
| `plots` | What does one run look like? (approximate model) |
| `compare_physical` | Does the tuning survive the *real* kinematics? **Run after every retune.** |
| `compare_gain` | Fixed vs. adaptive gain; how does the admissible gain depend on epsilon? |
| `compare_V` | Which Lyapunov function converges fastest within the velocity limits? |
| `search_k(kmax)` | Which frequency multipliers (k12, k112, k1112) allow the largest gain? |
| `convergence_map` | For which (epsilon, gamma) does the controller converge at all? |
| `derive_symbolic` | The algebra: Lie brackets, det G = 1, amplitudes a(z) (needs Symbolic Toolbox) |
| `trailers2` | Grushkovskaya's own log-V script (17 Oct 2026), kept as she sent it |

## Layout

```
2 trailers/
  params.m               every number of the study, with the tuning notes
  setup_paths.m          puts lib/ on the path (params.m calls it for you)
  selftest.m             regression test with stored reference values
  plots.m  compare_*.m  search_k.m  convergence_map.m  derive_symbolic.m
  trailers2.m            Grushkovskaya's script (self-contained copy of the controller)
  lib/
    simulate.m           THE closed-loop engine (approx | exact, fixed | adaptive gain)
    ctrl_frame.m         the numbers the controller freezes at the start of an interval
    u_ctrl.m             u1(t), u2(t) inside the interval
    amps.m               the five amplitudes a(z) = -h G^-1 grad Q
    control_bound.m      worst-case |u1|, |u2|
    gamma_max.m          largest gain that respects the velocity limits
    fields_approx.m      vector fields of the nilpotent approximation   (g1, g2)
    fields_exact.m       vector fields of the real kinematics            (f1, f2)
    lyapunov_h.m         quadratic / log / sqrt Lyapunov functions
    offhooked_coeffs.m   alpha, kappa, theta, beta and the frame M0 from the geometry
    x_to_z.m  z_to_x.m   physical state <-> privileged coordinates
    check_transform.m    self-check of that transform
    set_epsilon.m        change epsilon AND omega together
    satisfies_c1.m       resonance condition C1 of the theorem
    print_setup.m        one-screen summary of the current parameters
  _old/  _derivation/    retired files and the slow symbolic derivation of x -> z
```

## How the pieces fit

```
params.m --> P (struct) --+--> simulate(P, model, gain) --> out.t, out.z, out.u1, out.u2 ...
                          |        per sampling interval:
                          |          x --x_to_z--> z  (exact model only)
                          |          ctrl_frame(z)  : amplitudes frozen once
                          |          ode45 with u_ctrl(t) and fields_approx / fields_exact
                          +--> control_bound, gamma_max   (no simulation needed)
```

* **`u1`, `u2` are physical** (driving velocity in m/s, steering velocity in rad/s)
  and identical in both descriptions. Only the *state* is converted (x -> z).
  On a TurtleBot3: `u1 -> cmd_vel.linear.x`, `u2 -> cmd_vel.angular.z`.
* **Approximate vs exact.** The controller is designed on the nilpotent
  approximation `z' = u1 g1(z) + u2 g2(z)`; the robot obeys `x' = u1 f1(x) + u2 f2(x)`.
  `simulate(P)` runs the former, `simulate(P, 'exact')` the latter, with the same controller.
* **Angles** (`P.x0`) are in **radians**.

## Rules that save you from silent mistakes

1. **Change epsilon with `P = set_epsilon(P, eps)`**, never `P.epsilon = ...` (omega would go stale).
2. **Change the Lyapunov function in `params.m`, section 4** (`lyapunov_h(...)`). Nothing else has to move.
3. **Validate every retune with `compare_physical`.** A large epsilon looks fast on the
   approximate model and can diverge on the exact one (see the notes in `params.m`).
4. **Run `selftest` after every edit.** If you change parameters on purpose it will fail on
   purpose -- check the new numbers, then update the reference table at its top.
5. Never add `_old/` to the MATLAB path -- it holds retired files with the same names as the ones in `lib/` (old signatures).
6. Start states have to lie inside the region where the approximation is valid
   (for one trailer measured as roughly 8 cm / 18 deg, see the thesis notes `OBALKA_1PRIVES.md`); outside it a trailer winds up.

## Old file names -> new

| before (17 Oct 2026) | now |
|---|---|
| `calculations.m` (symbolic part) | `derive_symbolic.m` |
| `calculations.m` (simulation part) + `plots.m` | `plots.m` |
| `simulate_case.m`, `simulate_case_adaptive.m`, `simulate_case_physical.m` | `lib/simulate.m` (`'approx'`/`'exact'`, `'fixed'`/`'adaptive'`) |
| `u_ctrl.m`, `amps.m`, `control_bound.m`, `gamma_max.m` | `lib/` (same formulas; now take the struct `P`) |
| `g1_fun`, `g2_fun` pasted into four scripts | `lib/fields_approx.m` |
| `rhs_x` inside `simulate_case_physical.m` | `lib/fields_exact.m` |
| `satisfies_c1` inside `search_k.m` | `lib/satisfies_c1.m` |
| `P.delta`, `P.h_fun = @(z) ...` | `[P.h_fun, P.lyap_name] = lyapunov_h('sqrt', 0.005)` |

Everything else (`params.m`, `plots.m`, `compare_physical.m`, `compare_gain.m`,
`compare_V.m`, `convergence_map.m`, `search_k.m`) keeps its name.
The previous versions of all files are in `_old/before_refactor_2026-10-17/`.

## Verification

The refactor was checked by running the old and the new code side by side in GNU Octave
and comparing the numbers (coefficients, control bounds, gamma_max, and 12 s closed-loop
runs of all three variants agree to 12 digits; `search_k(15)` returns the identical table).
Plotting commands and `derive_symbolic.m` could not be run there; the latter was
re-checked in Python/SymPy (det G = 1, a(z) agrees with `lib/amps.m` to 1e-22).
`selftest` repeats the numerical part on your MATLAB.
