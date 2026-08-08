function [A_loc, b_loc, lo, hi] = ieee13_local_capability(community, up_frac, q_frac)
%IEEE13_LOCAL_CAPABILITY Local (per-inverter) box capability constraints
%A_loc*x_F <= b_loc, in the SAME 14-column order as A_F_sh/derive_shared_
%constraint.m: per player [645,611,652,671], each as [p_phases, q_phases].
%
%   [A_loc, b_loc, lo, hi] = IEEE13_LOCAL_CAPABILITY(community, up_frac, q_frac)
%
%   EXPLICIT, TRANSPARENT ASSUMPTION (no capability data was specified --
%   this is a first, clearly-labeled modeling choice, open to revision):
%
%   Each inverter is dispatched in S'' at its nameplate real power P_DER,i
%   (as fixed in ieee13_der.m) and is assumed MODESTLY OVERSIZED on the AC
%   side (a common PV "inverter loading ratio" > 1 in practice), giving it
%   some headroom for BOTH a real-power boost and reactive support:
%
%       Delta P_i (net consumption) in [-up_frac * P_DER_i,  +P_DER_i]
%           (can boost output by up_frac above nameplate; can fully curtail
%            down to zero output, i.e. increase net consumption by P_DER_i)
%       Delta Q_i               in [-q_frac * P_DER_i, +q_frac * P_DER_i]
%           (symmetric reactive support/absorption capability)
%
%   Defaults: up_frac = 0.20 (20% real-power headroom), q_frac = 0.30 (30%
%   reactive capability relative to nameplate real power). These are BOX
%   (independent P/Q) bounds, an approximation to the true circular/
%   elliptical S_rated capability curve of a real inverter -- adequate for
%   a first linear feasibility pass, not a substitute for a manufacturer
%   capability curve.
%
%   OUTPUTS
%   A_loc, b_loc : 28 x 14 / 28 x 1, such that A_loc*x_F <= b_loc encodes
%                  lo <= x_F <= hi.
%   lo, hi       : 14 x 1 bounds themselves, for reporting.

    if nargin < 2 || isempty(up_frac); up_frac = 0.20; end
    if nargin < 3 || isempty(q_frac);  q_frac  = 0.30; end

    Sbase = 5e6;
    P_DER_kW = [community.p_kW];          % 7-vector, community column order
    P_DER_pu = P_DER_kW(:) * 1e3 / Sbase;

    lo_p = -up_frac * P_DER_pu;
    hi_p =  P_DER_pu;
    lo_q = -q_frac * P_DER_pu;
    hi_q =  q_frac * P_DER_pu;

    lo = [lo_p; lo_q];
    hi = [hi_p; hi_q];

    nF = numel(lo);
    A_loc = [eye(nF); -eye(nF)];
    b_loc = [hi; -lo];
end
