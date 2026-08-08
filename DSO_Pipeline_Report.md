# A DSO-Side Physical Pipeline for the IEEE-13 Feeder: Construction, Validation, and the Local Voltage-Security Interface

## Abstract

We construct a clean, three-stage physical pipeline on the IEEE-13 three-phase
unbalanced test feeder,
$$
\mathcal S \;\rightarrow\; \mathcal S' \;\rightarrow\; \mathcal S'',
$$
progressing from the original benchmark to an active distribution system with
inverter-based DER at four future energy-community (EC) nodes, to an active
and congested system under a $\lambda=1.8$ load-stress factor. At every stage
we solve the Bolognani–Dörfler (1ACPF) tangent-plane manifold linearization
and independently validate against a custom Newton–Raphson nonlinear
three-phase AC solver. From the fully validated $\mathcal S''$ operating
point we derive the local input–output sensitivity map, the grid-wide
voltage-security constraint, and its projection onto the four EC
participants, $A_{F,sh}\Delta s_F \le b_{F,sh}$. We then run a pure
feasibility/minimum-effort analysis of that constraint (no game-theoretic
objective yet) and re-validate the resulting corrective dispatch in the
nonlinear AC model, and separately quantify how the linearization error of
the interface itself grows with distance from the $\mathcal S''$ anchor.
One methodological finding is reported in detail because it changed a
concrete result: the DSO interface must be anchored at the nonlinear
AC-converged operating point, not at the fast one-shot manifold estimate —
an anchoring bug that silently corrupted the exported constraint's RHS by
up to 0.19 p.u. was found and fixed during this work, verified by a
hard invariant check now built into the pipeline.

## 1. Motivation

The end goal (outside the scope of this report) is a generalized Nash
equilibrium game among four flexible energy-community participants,
constrained by the distribution network's voltage-security limits. Before
any game-theoretic formulation is introduced, the DSO side of that
interface needs to be established on its own: a physically meaningful,
independently validated operating point, from which a *local* network
constraint is derived and handed to the community layer. This report covers
exactly that DSO-side construction — nothing about player objectives, cost
structure, or equilibrium selection is addressed here.

## 2. Method

### 2.1 The Bolognani–Dörfler / 1ACPF manifold

The repository's existing machinery (`Nmatrix.m`, `Rmatrix.m`, `bracket.m`,
`rw.m`) implements the implicit linearization of the power-flow manifold
from Bolognani & Dörfler (Allerton 2015). Two forms of this linearization
are used, for two different purposes:

- **One-shot / flat-anchored** (`bd_linear_solve.m`): the original paper's
  method — a *single* linear solve, without iteration, using the flat/
  no-load reference voltage as the linearization point. This is the "fast"
  manifold solve; we use it at every stage purely to obtain a quick
  approximate operating point $x^\star$ and to measure how much error that
  approximation carries relative to the true nonlinear solution.
- **General / arbitrary-point** (`bd_tangent_general.m`): the same
  underlying tangent-plane construction, generalized to linearize around
  *any* point $(v^\star,\theta^\star)$ on or near the manifold:
$$
A_{x^\star} = \big[\,M(u^\star)R(u^\star),\; I\,\big],\qquad
M(u^\star)=\langle \mathrm{diag}(\overline{Yu^\star})\rangle + \langle
\mathrm{diag}(u^\star)\rangle N\langle Y\rangle,
$$
  where $\langle\cdot\rangle$ is the real-imaginary bracket operator, $N$
  the sign-flip operator on the $[v;\theta]$/$[p;q]$ split, and
  $R(v,\theta)$ the polar-to-rectangular rotation. This is the form that
  produces the final DSO-to-community interface, evaluated at the nonlinear
  AC operating point of $\mathcal S''$ (Section 5).

### 2.2 Nonlinear AC ground truth

No solver for unbalanced three-phase nonlinear AC power flow existed in
this repository (MATPOWER, used elsewhere in the codebase, is
positive-sequence/balanced only). `solve_acpf_3ph_rect.m` is a
purpose-built Newton–Raphson solver in rectangular coordinates, operating
directly on the same multiphase admittance matrix $Y$ built by `ieee13.m`,
used identically at every stage as the independent validation reference.
It converges in 3–4 Newton iterations to a mismatch below $10^{-12}$ p.u.
in every case reported here.

### 2.3 Physical connectivity mask

The original feeder mixes fully three-phase and partially connected
buses. Because `ieee13.m`'s admittance matrix regularizes missing phases
with a nonzero dummy self-impedance (needed to keep $Y$ invertible even on
laterals with no physical conductor on some phase), the true missing-phase
structure is *not* recoverable from $Y$ alone. It is instead taken from the
NaN pattern of the benchmark's own published solution table (Kersting), a
single source of truth reused at every stage since topology and phase
connectivity never change across $\mathcal S,\mathcal S',\mathcal S''$: 29
of 36 non-slack bus-phase terminals are physically connected.

## 3. Stage definitions

### 3.1 $\mathcal S$ — original benchmark (unmodified)

Exactly as provided: same topology, loads, capacitor banks, slack
conditions, and voltage limits ($0.95\le|V|\le1.05$ p.u.) as the published
IEEE-13 feeder. No modification.

### 3.2 $\mathcal S'$ — active distribution system (+DER)

Inverter-based DER is added at the four future EC nodes — **645, 611, 652,
671** — using the additive, minimal-modification convention
$$
P_i^{\text{net}} = P_i^{\text{load}} - P_i^{\text{DER}}, \qquad
Q_i^{\text{net}} = Q_i^{\text{load}} \;\;(\,Q^{\text{DER}}=0\,),
$$
i.e. unity-power-factor inverters placed *alongside* the existing load
rather than replacing it, so no other part of the feeder is touched. DER is
placed only on bus-phases already physically connected *and* already
carrying load (645: phases b, c; 611: phase c; 652: phase a; 671: phases
a, b, c), sized at roughly 1.5–2$\times$ the local phase load:

| bus | phase(s) | DER [kW] |
|---|---|---|
| 645 | b, c | 300, 150 |
| 611 | c | 300 |
| 652 | a | 250 |
| 671 | a, b, c | 500, 500, 500 |

Total DER: **2500 kW**, against an original total feeder load of
$\approx$3466 kW — a substantial DG penetration. These sizes and the
additive convention are explicit, documented modeling choices (there is no
external capacity spec to match against), open to revision.

### 3.3 $\mathcal S''$ — active and congested system

A $\lambda=1.8$ load-side stress is applied **starting from $\mathcal
S'$**, scaling only the demand component of net consumption:
$$
S_{\text{load}}'' = \lambda\, S_{\text{load}}(\mathcal S), \qquad
S_{\text{cap}}'' = S_{\text{cap}}(\mathcal S)\ \text{(fixed)}, \qquad
S_{\text{DER}}'' = S_{\text{DER}}(\mathcal S')\ \text{(fixed, not scaled)}.
$$
Capacitor banks and the $\mathcal S'$ DER dispatch are explicitly *not*
scaled by $\lambda$ — per the pipeline's stress convention, an increase in
feeder loading is a demand-side phenomenon; scaling fixed devices or
already-installed generation by the same factor would not represent that
physical scenario. Total feeder active load rises from 3466 kW to 6239 kW;
DER stays fixed at 2500 kW.

## 4. Per-stage results

At every stage: (1) one-shot manifold solve, (2) fresh tangent-plane
linearization at *that stage's own* operating point (never reused from a
previous stage), (3) independent nonlinear AC solve, (4) comparison.

| stage | $V_{min}$ [pu] | $V_{max}$ [pu] | non-slack violations (under/over) | one-shot lin. vs AC: $e_{max}$ | $e_{RMS}$ |
|---|---|---|---|---|---|
| $\mathcal S$   | 0.9737 (611.c) | 1.0557 | 0 / 4 | 5.64e-3 | 3.93e-3 |
| $\mathcal S'$  | 1.0005 (675.a) | 1.0603 | 0 / 4 | 2.17e-3 | 1.14e-3 |
| $\mathcal S''$ | 0.9051 (675.c) | 1.0607 | 14 / 4 | 1.93e-2 | 1.15e-2 |

$\mathcal S$'s nonlinear AC solution matches the published Kersting
benchmark to $3.8\times10^{-3}$ p.u. max error, confirming the base model
and solver are correct before anything is modified.

![Voltage profile: S vs S' vs S'', nonlinear AC, all three phases](figures/fig1_voltage_profile.png)

![Violation counts per stage](figures/fig2_violations.png)

![One-shot linearization error per stage](figures/fig3_linearization_error.png)

**Reading these results:**

- **$\mathcal S$ is weakly, not trivially, constrained.** Four non-slack
  overvoltage violations already exist (671.b, 692.b, 675.b, 680.b, all
  $\approx$1.053–1.056 p.u.), caused by the fixed 675 capacitor bank
  combined with structurally lighter phase-b loading in this unbalanced
  feeder — not something the later stages introduce.
- **DER in $\mathcal S'$ raises the whole profile.** $V_{min}$ rises from
  0.9737 to 1.0005 p.u. (611, the feeder's weakest point, benefits most
  from reduced net consumption upstream of it), but the same four
  overvoltage violations persist, nudged slightly *higher* — DER has no
  mechanism to relieve a phase-b-specific overvoltage caused by loading it
  never touches on phases a/c.
- **$\mathcal S''$ is heavily stressed, with opposing pressures.**
  18 of 29 non-slack terminals violate: undervoltage on phases a/c across
  most of the feeder, while the phase-b overvoltage persists unchanged.
  This is a materially richer congestion pattern than a single-direction
  "everything is low" scenario — a corrective action cannot simply "raise
  voltage everywhere."
- **The one-shot linear approximation degrades under stress**, from 0.56%
  max error at $\mathcal S$ to 1.93% at $\mathcal S''$ — and the error is
  not merely cosmetic: at $\mathcal S''$ the one-shot model reports only 1
  non-slack overvoltage violation where the true AC solution has 4,
  changing the inferred constraint set's row count, not just its precision.

## 5. The DSO-to-community interface

### 5.1 Anchoring — and a bug found in the process

Section 4's derivation (`derive_shared_constraint.m`) partitions the
tangent-plane system into slack/non-slack blocks,
$$
A_z \Delta z_{ns} + A_s \Delta s_{ns} = 0 \;\Rightarrow\;
\Delta z_{ns} = S_{zs}\Delta s_{ns} = -A_z^{-1}A_s\,\Delta s_{ns}, \qquad
S_{zs}=\begin{bmatrix}S_v\\S_\theta\end{bmatrix},
$$
then builds the grid-wide constraint at the *selected* operating point
$v^{\star\prime\prime}$:
$$
A_{grid}=\begin{bmatrix}S_v\\-S_v\end{bmatrix}, \qquad
b_{grid}=\begin{bmatrix}v_{max}-v^{\star\prime\prime}\\v^{\star\prime\prime}-v_{min}\end{bmatrix}.
$$

**An initial implementation anchored $v^{\star\prime\prime}$ at the
one-shot manifold estimate rather than the nonlinear AC-converged point.**
Because the pipeline explicitly adopts $\mathcal S''$ "as the operating
system" only *after* the manifold-vs-AC comparison (Section 4), the
physically selected $v^{\star\prime\prime}$ is the validated AC solution —
anchoring at the faster but less accurate estimate instead silently changed
which rows the exported constraint reported as violated (1 upper / 14 lower
instead of the correct 4 upper / 14 lower). A second, independent bug
compounded this: the projection step included a spurious RHS term
$-A_F\bar s_F-A_P\bar s_P$ inherited from an earlier draft, when the
$\Delta s_P=0$ deviation formulation requires $b_{F,sh}=b_{grid}$ exactly.
Together these corrupted the exported margins by up to 0.19 p.u. Both are
fixed, and a hard invariant assertion is now built into the pipeline: it
refuses to export unless the constraint's violation pattern exactly
reproduces the independently computed nonlinear AC ground truth. This is
recorded here because it is a natural mistake for this kind of derivation
to make silently, and the general lesson — *derive the constraint from the
same point you validated, and check it against ground truth before
trusting it* — is broadly applicable.

### 5.2 The resulting interface

With the corrected anchor at $\mathcal S''$'s nonlinear AC solution:

$$
A_{grid}\in\mathbb R^{58\times72}\ \ (29\text{ physical non-slack terminals}), \qquad
A_{F,sh}\in\mathbb R^{58\times14}, \quad b_{F,sh}\in\mathbb R^{58},
$$

projected onto the 14 controllable $[p,q]$ variables of the four
participants (645, 611, 652, 671, in that order). Minimum voltage margins:
upper $-0.0107$ p.u., lower $-0.0449$ p.u. — 4 upper-margin and 14
lower-margin rows are negative, exactly matching the AC ground truth's 4
overvoltage / 14 undervoltage count. Exported to
`data/network_to_dynect_congested.mat`.

![|S_v| sensitivity heatmap: community [p,q] columns vs physical terminal rows](figures/fig4_sensitivity_heatmap.png)

The sensitivity heatmap shows the expected structure: the strongest
couplings are concentrated among electrically close terminals (671, 692,
675, 680, 684, 611, 652 — all downstream of the 671 lateral where three of
the four community nodes sit), with 645's phases coupling more broadly
across the 632-side of the feeder. $q645b$ and $p671b/p671c$ stand out as
particularly influential columns.

## 6. Pure feasibility check (no game objectives)

Before any game-theoretic formulation, we test only
$$
\mathcal F=\{x_F:\ A_{F,sh}x_F\le b_{F,sh},\ A_{loc}x_F\le b_{loc}\}.
$$
$A_{loc}$ encodes explicit, documented per-inverter capability bounds (no
manufacturer data was available; this is a first, revisable assumption):
real power boostable by 20% above nameplate and fully curtailable, reactive
support/absorption up to 30% of nameplate (`ieee13_local_capability.m`).

- **$x_F=0\notin\mathcal F$**, confirmed directly from $b_{F,sh}$'s sign
  pattern (4 upper / 14 lower negative rows), matching the AC ground truth.
- **$\mathcal F\neq\emptyset$**: a feasibility LP (`glpk`) finds a strictly
  feasible point, $\|x_{feas}\|_2=0.0924$ p.u.
- **Minimum-effort dispatch** $x_{min}=\arg\min\tfrac12\|x_F\|_2^2$ s.t.
  $\mathcal F$ (solved via ADMM — no QP solver package was available in
  this environment; convergence and constraint satisfaction were verified
  independently rather than trusted from the solver's own residuals):
  $\|x_{min}\|_2=0.0552$ p.u., resolving all 18 rows *in the linear
  tangent model*. Several local capability bounds are saturated at this
  optimum (q645b, p652a, p671b, p671c, q671c all sit exactly at their
  assumed limits) — the assumed capability set is not comfortably interior
  to the optimum, it is close to binding.
- **Nonlinear AC re-validation of $x_{min}$: not fully certified.**
  $V_{min}$ improves from 0.9051 to 0.9234 p.u. and all 4 overvoltage
  violations clear completely, but 9 of the 14 undervoltage violations
  remain (634.a/c, 671.c, 692.c, 675.a/c, 680.c, 684.c, 611.c).

![S'' baseline vs minimum-effort-corrected voltage profile (nonlinear AC)](figures/fig5_correction_before_after.png)

This does **not** invalidate the derived interface — it is a statement
about the *range of validity* of a local (tangent-plane) constraint under a
large correction, which Section 7 quantifies directly.

## 7. Quantifying the input–output map's accuracy away from the anchor

The affine map $v_{ns}\approx v^{\star\prime\prime}+S_v\Delta s_{ns}$ is
exact at $\Delta s=0$ by construction; how it degrades away from that point
had not been measured directly. `run_io_map_accuracy.m` sweeps 8 random
directions in the 14-dimensional community flexibility space at 6
magnitude fractions (5%–100%) of the assumed local capability box, plus the
$x_{feas}$/$x_{min}$ points from Section 6, comparing the tangent-map
prediction against a **fresh, independent** nonlinear AC solve at each
perturbed point (not the AC solve used to build the tangent plane).

![IO-map error vs distance from the S'' anchor (log-log)](figures/fig6_io_map_accuracy.png)

| point | $\|\Delta s_F\|_2$ [pu] | $e_{max}$ [pu] | $e_{RMS}$ [pu] |
|---|---|---|---|
| anchor | 0 | $1.1\times10^{-16}$ (verified, not assumed) | 0 |
| smallest sweep sample | 0.00053 | $5.3\times10^{-8}$ | $2.3\times10^{-8}$ |
| largest sweep sample (100% of box) | 0.0482 | $6.8\times10^{-4}$ | $3.7\times10^{-4}$ |
| $x_{min}$ | 0.0552 | $9.1\times10^{-4}$ | $4.2\times10^{-4}$ |
| $x_{feas}$ | 0.0924 | $5.3\times10^{-3}$ | $2.6\times10^{-3}$ |

The error grows approximately quadratically with distance from the anchor
— consistent with a second-order Taylor remainder, now empirically
confirmed rather than assumed. Within the assumed capability box under
generic random perturbations, error stays small ($\le7\times10^{-4}$
p.u.). $x_{feas}$ and $x_{min}$ are corner-like points of the box (several
capability bounds simultaneously active), which is exactly why their norms
and errors exceed what a single random direction reaches within the same
box: this is the direct, now-quantified explanation for Section 6's partial
AC-certification result — not a separate anomaly.

## 8. Discussion and open items

**Established by this work:**
1. A validated, three-stage physical progression $\mathcal S\to\mathcal
   S'\to\mathcal S''$, each independently checked against a purpose-built
   nonlinear three-phase AC solver.
2. A correctly anchored, invariant-checked DSO-to-community interface
   $A_{F,sh},b_{F,sh}$, derived only from $\mathcal S''$.
3. Confirmation that the community's local constraint is non-trivial
   ($x_F=0\notin\mathcal F$) and, under the assumed capability set,
   feasible ($\mathcal F\neq\emptyset$) in the linear model.
4. A quantified, not assumed, characterization of the interface's local
   accuracy: exact at the anchor, small within the assumed capability
   region, and measurably degraded at the specific corner-like points a
   real optimization would select.

**Not yet established, and should not be assumed:**
- That the four EC participants, under the current capability assumption,
  can *fully* restore nonlinear voltage security at $\lambda=1.8$ — Section
  6 shows a partial (9/18) result for one dispatch, not a general claim.
- Whether a materially different, but still transparent, local capability
  assumption (e.g. larger reactive headroom) changes this outcome —
  untested.
- How the fixed 2500 kW of DER, if redistributed differently across the
  four nodes (still summing to 2500 kW), reshapes $S_v$ and hence the
  later game's structure — the natural next study, deliberately deferred
  until this checkpoint passed.

**Two natural next steps**, not yet decided between:
1. **SLP-style iteration**: apply a correction, re-solve AC, re-linearize
   at the new point, re-solve feasibility/min-effort there, repeat.
2. **Revisit the physical scenario**: loosen the (explicitly documented,
   unmeasured) capability assumption, or reduce $\lambda$, and re-run the
   same checkpoint.

## 9. Reproducibility

All results in this report are regenerated by, in order:
```
octave --no-gui --eval "run_ieee13_pipeline"       % Sections 3-5
octave --no-gui --eval "run_feasibility_check"     % Section 6
octave --no-gui --eval "run_io_map_accuracy"       % Section 7
octave --no-gui --eval "generate_report_figures"   % Figures 1-5 (fig6 is copied from run_io_map_accuracy's output)
```
Each script re-derives everything it needs from `ieee13.m` onward rather
than depending on a previous run's saved state, so any of them can be run
standalone after a `git pull`. Verified in Octave 8.4 and MATLAB (both
independently, same numerical results to solver-precision).

## Appendix: file map

| file | role |
|---|---|
| `ieee13.m`, `Nmatrix.m`, `Rmatrix.m`, `bracket.m`, `rw.m` | original paper's machinery, unmodified |
| `solve_acpf_3ph_rect.m` | nonlinear 3-phase AC solver (new, this work) |
| `ieee13_der.m`, `ieee13_congest.m`, `ieee13_capacitors.m` | $\mathcal S\to\mathcal S'\to\mathcal S''$ construction |
| `bd_linear_solve.m`, `bd_tangent_general.m` | one-shot and general tangent-plane linearizations |
| `voltage_report.m`, `compare_voltages.m` | shared diagnostics used identically at every stage |
| `run_stage.m` | per-stage orchestration |
| `derive_shared_constraint.m` | Sections 4-6 of the pipeline spec: $S_v$, $A_{grid}$, $A_{F,sh}$/$b_{F,sh}$ |
| `run_ieee13_pipeline.m` | main driver: $\mathcal S,\mathcal S',\mathcal S''$ + interface export |
| `ieee13_local_capability.m`, `feasibility_lp.m`, `min_effort_qp.m`, `run_feasibility_check.m` | Section 6 |
| `run_io_map_accuracy.m` | Section 7 |
| `generate_report_figures.m` | Figures 1-5 |
| `legacy/` | superseded exploratory version, kept for reference |
