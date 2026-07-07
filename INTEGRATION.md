# GNE input-output map layer — how it sits on top of `ETH_Manifold_1ACPF`

This is a **derivation layer**, not a fork. It reuses the repo's tangent-manifold
machinery unchanged and post-processes it into an explicit affine map
`y = S·u + m` and linear constraints `A_sh·x ≤ b_sh`. It implements the bridge
document (`context_ACPF_to_shared_constraint_bridge.md`), §5–§14.

## What is reused, untouched
`bracket.m`, `Rmatrix.m`, `Nmatrix.m`, `rw.m` — the exact building blocks the
paper's code uses. No edits.

## The one insight
`example_case14.m` already assembles

```
Agrid = [(JJ + UU*NN*YY)*PP  -eye(2n)]   ==   A_{ξ*} = [A_vd  -I]
```

It then stacks bus-model rows `AA/BB` that **pin the injections to numbers** and
solves once (the "1ACPF" approximate operating point). This layer instead keeps
the injections **symbolic**: it partitions `A_vd` into unknown-voltage columns
(`A_U`) vs input-power rows and inverts, giving `S`. So §5–§8 of the bridge is
literally "drop the `AA/BB` block and do the partition."

## One correction baked in
`example_case14.m` hardcodes `PP = Rmatrix(ones,zeros)` because it linearizes at a
**flat** voltage. This layer linearizes at the **true AC operating point** `ξ*`
(bridge §14 step 2), so it uses `PP = Rmatrix(V0,A0)`. At a non-flat point the
flat form injects ~1e-3 error into every entry of `S` (verified). See `build_Avd.m`.

## Files added
| File | Bridge § | Role |
|---|---|---|
| `build_Avd.m` | 2 | tangent matrix `A_vd`, factored out + generalized to any `ξ*` |
| `build_io_map.m` | 5–7, 9 | partition → `S`, `m` (PQ-for-everyone v1 scope) |
| `build_shared_constraint.m` | 8, 10, 11 | `A_sh_grid`, `b_sh'`, per-player `Π_i` blocks |
| `nrpf.m` | (test util) | standalone Newton-Raphson AC PF, no MatPower |
| `validate_io_map.m` | 12 | Tests 1–4, run **before** the GNE solver |
| `example_gne_4bus.m` | 14 | end-to-end self-contained demo |

## Run it
```
octave --no-gui example_gne_4bus.m
```
Self-contained (manual Ybus, no MatPower). Expected: all four tests PASS, with
Test 2 (finite-difference vs true AC) at ~1e-10.

## Two v1-scope guardrails from the bridge (§9, §13) enforced in code
- **PQ-for-everyone.** Every non-slack bus enters as a controllable PQ injection;
  no PV buses in the tangent map. `build_io_map` builds `A_U` from the non-slack
  `[V;θ]` columns only. "Game technology PV" ≠ "ACPF PV bus" — never conflate.
- **γ never enters physics.** `Π_i = [[1,0,0];[0,1,0]]` drops `γ_i`; Test 4 asserts
  the third column of every player block is identically zero.

## Wiring to MatPower / `case33bw` (the real target)
Replace the manual Ybus block in the driver with:
```
mpc = loadcase('case33bw'); results = runpf(mpc, mpoption('OUT_ALL',0));
Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
V0 = results.bus(:,VM); A0 = results.bus(:,VA)*pi/180;   % TRUE operating point
```
Then `build_io_map(Ybus, V0, A0, slack, mon_buses)` unchanged. `nrpf.m` is only
for the self-contained test; with MatPower present, Test 2 can call `runpf`
instead for an even stronger check.

## Hand-off
`sh.A_sh_blocks` (one `m_sh × 3` matrix per FL player) and `sh.b_sh` are exactly
what `StaticGNEGame` needs (bridge §11). Final constraint:
`Σ_i A_sh_blocks{i}·x_i ≤ b_sh`.
