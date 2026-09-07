function [p2, q2, meta] = ieee13_congest(p_S, q_S, q_cap, p_der, q_der, lambda)
%IEEE13_CONGEST Apply a load-side congestion/stress factor to the active
%system S', producing the stressed active system S''.
%
%   [p2, q2, meta] = IEEE13_CONGEST(p_S, q_S, q_cap, p_der, q_der, lambda)
%
%   CONVENTION (Section 3 of the DSO-pipeline spec)
%   ------------------------------------------------------------
%   The stress represents increased feeder loading, so the physically
%   appropriate transformation scales the LOAD component of demand only:
%
%       S_load'' = lambda * S_load        (S_load taken from the ORIGINAL
%                                           benchmark S -- p_S, q_S+q_cap --
%                                           since S' only added DER on top
%                                           of that same load, it never
%                                           changed it)
%
%   Fixed devices are preserved unless there is a physical reason to change
%   them:
%     * capacitor banks (q_cap) are NOT scaled;
%     * the DER operating points introduced in S' (p_der, q_der) are NOT
%       automatically scaled by lambda -- they are a separate, deliberately
%       chosen physical quantity (existing inverter capacity does not grow
%       just because neighboring load grows). They are simply re-subtracted
%       from the newly scaled load to get the S'' net consumption:
%
%       P_i^net,S''  = lambda * P_i^load,S   - P_i^DER
%       Q_i^net,S''  = lambda * (Q_i^load,S + Q_cap_i) - Q_cap_i - Q_i^DER
%
%   INPUTS
%   ------
%   p_S, q_S     : 3n x 1 original-benchmark net consumption (from ieee13()).
%                  p_S is pure load (no capacitor component); q_S already
%                  has the capacitor banks netted in (q_S = q_load - q_cap).
%   q_cap        : 3n x 1 fixed capacitor injections (ieee13_capacitors.m).
%   p_der, q_der : 3n x 1 DER capacity introduced in S' (ieee13_der.m).
%   lambda       : scalar load-stress factor (e.g. 1.8).
%
%   OUTPUTS
%   -------
%   p2, q2 : 3n x 1 net consumption of S''.
%   meta   : bookkeeping struct (lambda, load totals before/after).

    q_load_S = q_S + q_cap;

    p2 = lambda * p_S - p_der;
    q2 = lambda * q_load_S - q_cap - q_der;

    % Slack rows (1:3) are NaN in ieee13()'s p,q (no load spec there);
    % totals must be summed over non-slack rows only.
    nonSlack = 4:numel(p_S);

    meta.lambda = lambda;
    meta.total_P_load_S  = sum(p_S(nonSlack));
    meta.total_P_load_S2 = lambda * sum(p_S(nonSlack));
    meta.total_Q_load_S  = sum(q_load_S(nonSlack));
    meta.total_Q_load_S2 = lambda * sum(q_load_S(nonSlack));
    meta.convention = ['S_load'''' = lambda*S_load(S); capacitors fixed; ', ...
        'DER (from S'') held fixed, not scaled by lambda'];
end
