1ACPF
=========

This repository contains the Matlab / GNU Octave and MatPower source code used in the paper:

S. Bolognani and F. Dörfler, "Fast power system analysis via implicit linearization of the power flow manifold,"
presented at the 53rd Annual Allerton Conference on Communication, Control, and Computing in Monticello, IL (2015).

- MatPower is available at http://www.pserc.cornell.edu/matpower/
- GNU Octave is available at http://www.gnu.org/software/octave/

## DSO-side physical pipeline (IEEE-13, three-phase unbalanced)

`run_ieee13_pipeline.m` is the entry point for the DSO-side pipeline built on
top of the original paper's method. It constructs three distinct physical
systems in sequence -- `S` (original benchmark) -> `S'` (active feeder, +DER
at buses 645/611/652/671) -> `S''` (S' under a lambda=1.8 load-congestion
stress) -- validating each against an independent nonlinear three-phase AC
solver before deriving, from `S''` only, the local voltage-sensitivity map
and the DSO-to-energy-community shared constraint `A_F_sh * Delta s_F <=
b_F_sh`. Run with:

```
octave --no-gui --eval "run_ieee13_pipeline"
```

See the header comments in `run_ieee13_pipeline.m` for the full stage
breakdown and `legacy/README.md` for notes on the earlier exploratory
version of this pipeline that it supersedes.

### DSO output (current result, lambda=1.8, nonlinear AC ground truth)

| system | Vmin [pu] | non-slack violations |
|---|---|---|
| S (original benchmark) | 0.9737 | 4 (all mild, phase-b overvoltage) |
| S' (+DER at 645/611/652/671) | 1.0005 | 4 (same 4, unchanged) |
| S'' (S' congested, lambda=1.8) | 0.9051 | 18 (14 under + 4 over) |

`S''`'s nonlinear AC-converged operating point is the adopted operating
system. From it, `derive_shared_constraint.m` produces the DSO-to-community
interface:

```
A_F_sh in R^(58x14),  b_F_sh in R^58
A_F_sh * Delta s_F <= b_F_sh
```

exported to `data/network_to_dynect_congested.mat`. This is a **local**
constraint -- valid as the tangent-plane linearization of the voltage-
security envelope around `S''`, which is its exact theoretical meaning; it
does not by itself claim anything about corrections far from that point.

`run_feasibility_check.m` is a **downstream diagnostic** on that interface,
not a precondition for it: it confirms `x_F=0` is infeasible (i.e. the
community constraint is non-trivial) and that the linear feasible set `F`
is non-empty, then separately checks whether one particular corrective
dispatch still holds up when re-solved in the true nonlinear AC model. A
negative answer there is a statement about the range of validity of the
local constraint under a large correction, not a retraction of `S''` or of
`A_F_sh`/`b_F_sh`.


