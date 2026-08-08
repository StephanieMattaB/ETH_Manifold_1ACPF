function [p_der, q_der, meta] = ieee13_der(n)
%IEEE13_DER Inverter-based DER injected at the four future energy-community
%nodes, defining the S -> S' modification of the DSO pipeline.
%
%   [p_der, q_der, meta] = IEEE13_DER(n)
%
%   CONVENTION (explicit, per the DSO-pipeline spec, Section 2)
%   ------------------------------------------------------------
%   ieee13.m's p,q are POSITIVE NET CONSUMPTION = P_load - P_generation.
%   DER is introduced ALONGSIDE the existing local demand (not by deleting
%   or replacing it), i.e.
%
%       P_i^net = P_i^load - P_i^DER          (consumption convention)
%       Q_i^net = Q_i^load - Q_i^DER
%
%   This is the minimal, physically interpretable modification: it adds a
%   behind-the-meter inverter at each of the four community buses without
%   reshaping any other part of the feeder (loads, capacitors, topology,
%   phase connectivity, and the slack bus are all left exactly as in S).
%
%   Reactive power is treated with the SAME convention and, for this first
%   active-system design, DER is modeled as unity-power-factor inverters:
%   Q_i^DER = 0, so Q_i^net = Q_i^load unchanged. (A non-zero Q_DER, e.g.
%   volt-var support, would be a later refinement -- not needed to obtain a
%   physically meaningful active/DG feeder.)
%
%   SIZING
%   ------
%   Only phases that are PHYSICALLY CONNECTED at each bus (per the original
%   missing-phase structure of S, preserved throughout) receive DER. Sizes
%   are chosen at roughly 1.5-2x the existing local phase load (or a
%   comparable round figure on phases with no local load), large enough to
%   create genuine reverse-flow / active-network behavior, without being
%   sized so large that S' is already congested (congestion is deliberately
%   introduced only at the next stage, S'').
%
%   bus 645 (idx 5):  phases b,c  ->  300 kW (b), 150 kW (c)
%   bus 611 (idx 12): phase c     ->  300 kW
%   bus 652 (idx 13): phase a     ->  250 kW
%   bus 671 (idx 7):  phases a,b,c -> 500 kW each
%
%   Total DER = 2500 kW (unity PF) against an original feeder total load of
%   ~3466 kW -- a substantial, but not implausible, DG penetration.
%
%   OUTPUT
%   ------
%   p_der, q_der : 3n x 1 vectors (per-unit on Sbase = 5e6 VA), zero
%                  everywhere except the community bus-phases above.
%                  Positive value = DER active/reactive CAPACITY at that
%                  bus-phase (to be SUBTRACTED from net consumption).
%   meta         : struct with the bus/phase/value bookkeeping, for
%                  reporting and for the S''-stage congestion step (which
%                  must NOT scale this vector by lambda).

    Sbase = 5e6;

    p_der = zeros(3*n, 1);
    q_der = zeros(3*n, 1);   % unity PF: no reactive DER capacity

    phase.a = 1; phase.b = 2; phase.c = 3;

    community = struct( ...
        'id',    {645,                611,           652,           671}, ...
        'bus',   {5,                  12,            13,            7}, ...
        'phase', {[phase.b phase.c],  phase.c,       phase.a,       [phase.a phase.b phase.c]}, ...
        'p_kW',  {[300 150],          300,           250,           [500 500 500]});

    for k = 1:numel(community)
        c = community(k);
        idx = rw(c.bus, c.phase);
        p_der(idx) = p_der(idx) + c.p_kW(:) * 1e3 / Sbase;
    end

    meta.community   = community;
    meta.player_ids  = [645; 611; 652; 671];
    meta.total_DER_kW = sum([community.p_kW]);
    meta.convention  = 'P_net = P_load - P_DER ; Q_net = Q_load (unity-PF inverters, Q_DER = 0)';
end
