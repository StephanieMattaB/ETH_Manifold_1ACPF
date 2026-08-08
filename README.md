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


