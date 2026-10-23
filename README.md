# MATLAB code of the bachelor thesis

Stabilization of a nonholonomic TurtleBot3 with off-hooked trailers by an
oscillating (Lie-bracket) controller on the nilpotent approximation
(Zuyev & Grushkovskaya, arXiv:2403.03341). Author: Markéta Hraničková.

## Where to start

1. Put **one** of the folders `1 trailer` / `2 trailers` on the MATLAB path, not both (each has its own
   `params.m`). Easiest: `cd` into the folder.
2. In `2 trailers/` (the main case) run  `trailer2('selftest')`  (about 30-45 s).
   It must end with `SELFTEST PASSED`; it compares the code with stored reference numbers.
3. `trailer2` alone lists all tasks; `trailer2('plots')` shows one closed-loop run.
4. In `1 trailer/` the same check is `trailer1('selftest')` (about 20 s).

All numbers (geometry, gain, period epsilon, k-multipliers, start, limits) are in
`params.m`. The only other place with numbers is the design block at the top of
`thesis_sim.m` / `thesis_figures.m`, see the table below.

## Folders and files

**Main code** (what the thesis is built on)

| path | what it is |
|---|---|
| `2 trailers/trailer2.m` | the whole two-trailer study in one file: simulators, controller, gains, vector fields. Call as `trailer2('task', ...)`. |
| `2 trailers/params.m` | every number of the two-trailer case (`P = params()` or `params(eps)`). Lyapunov function there: sqrt V (delta = 0.005). |
| `2 trailers/offhooked_coeffs.m`, `x_to_z.m` | geometry-only coefficients; physical state x -> privileged coordinates z. |
| `2 trailers/thesis_sim.m`, `thesis_figures.m` | data and figures quoted in the thesis (simulation chapter; feasibility limits). The design block at the top of each (eps, k, V) is written in the file itself, not read from `params.m`, so the chapter does not change when `params.m` is edited. |
| `1 trailer/` | the same for one trailer: `trailer1.m`, `params.m`, `thesis_sim.m`, `thesis_figures.m`, `derive.m` (symbolic, needs the Symbolic Math Toolbox). Lyapunov function everywhere in this folder: Victoria's log V, h = 1/(1 + c |z|^2), c = 10. |
| `thesis_figures_compare.m` | figure comparing r = 2, 3, 4; run after both `thesis_figures`. |
| `2 trailers/_derivation/x0_to_z0.m` | symbolic version of x -> z (slow, only for derivation). |

**Victoria's script**

| path | what it is |
|---|---|
| `2 trailers/trailers2.m` | Victoria's original script (log Lyapunov function), independent of `trailer2.m`. Kept unchanged for comparison. |

**Hyperexponential convergence** (h = |log ||z||^2|, task of 23 Oct 2026)

| path | what it is |
|---|---|
| `2 trailers/plots_hyper.m` | exact model, parameters of the paper, h = 1 against h = |log|; figures saved to `hyperexponential 2 trailers/`. |
| `1 trailer/plots_hyper_1tr.m` | one trailer, robot geometry, within the velocity limits: log V (c = 10) against h = |log|; figures saved to `hyperexponential 1 trailer/`. |

**Investigations** (`2 trailers/investigations/`; why the exact two-trailer model does not always converge)

Each file starts with a comment that says what it tests and how to call it. They add the folder `2 trailers` to the path themselves, so they can be run from anywhere.

| file | question |
|---|---|
| `compare_with_original.m` | does our exact-model simulator reproduce the older project `original_system` (not part of this repository)? |
| `bisect_regime.m` | walking from the settings of the old project to our design, which change first breaks convergence of the exact model? |
| `k_test.m` | does it converge for other multipliers k? |
| `angle_scan.m` | does the initial angle matter? |
| `paper_real.m` | the "real" version of the paper's simulation (paper parameters, exact kinematics). |

The cause of the non-convergence with the robot parameters is **not understood yet**; these scripts only narrow it down.

## Two models, one controller

* **approximation:** `z' = u1 g1(z) + u2 g2(z)` (nilpotent approximation), what the controller is designed on;
* **exact kinematics:** `x' = u1 f1(x) + u2 f2(x)` (`simulate_case_physical`), stands in for the real robot.

The same physical `u1, u2` drive both. At the start of every interval of length epsilon
the state x is converted to z, frozen, and the controls for that interval are computed from it.
Always check a design against the exact kinematics (`compare_physical`), not only on the approximation.

## How to change the Lyapunov function V

The controller sees V only through one scalar factor `h`, for a radial `V = phi(||z||)`:
`a(z) = -gamma * h(z) * G^-1(z) * z`,  `h(z) = phi'(r) / (gamma * r)`. So changing V means changing `h_fun`.

1. **Main place: `params.m` of the folder you work in.**
   * `1 trailer/params.m`: the lines `P.c_log`, `P.h_fun`, `P.lyap_name` (log V, c = 10). A sqrt alternative is in the comment below them.
   * `2 trailers/params.m`: the lines `P.delta`, `P.h_fun`, `P.lyap_name` (currently sqrt V, delta = 0.005).
2. **Retune `P.gamma`.** The admissible gain depends on V. For one trailer:
   `trailer1('gamma_ray', z0, alpha, omega, eps, k12, k112, 0.22, 2.84, h, 0.8)`;
   for two trailers `trailer2('compare_V')` prints the admissible gain of each candidate.
3. **`selftest` will fail on purpose** after the change, because its reference numbers belong to the old V: check the new
   numbers, then update the table in `selftest` (in `trailer1.m` / `trailer2.m`).
4. **Some files have their own copy of V and ignore `params.m`:**
   * `thesis_sim.m`, `thesis_figures.m` (both folders): the block `DES.*` / `h_use` at the top of the file,
   * the scripts in `2 trailers/investigations/` and `plots_hyper*.m` define their own `h`,
   * `2 trailers/trailers2.m` is the original script with its own log V.

## Things to know

* The optional phase offsets of the oscillating terms (`getappdata(0,'ctrl_phase')`, set only by
  `thesis_sim(..., struct('phase', [..]))`) are an **exploratory test outside the proof**. The default is zero,
  i.e. the controller of the paper.
* Hitch-angle limit: the printed "largest hitch angle" in `thesis_sim` should be compared with what the real joint
  allows (about 90 deg or less); the URDF value 1.40 rad is only a placeholder.
* `thesis_figures` and `thesis_sim` write PDFs and .mat files into `P.fig_dir`, which is the folder `figures/` next to
  `params.m` (i.e. `1 trailer/figures`, `2 trailers/figures`; created automatically). `thesis_figures_compare` reads from both.
