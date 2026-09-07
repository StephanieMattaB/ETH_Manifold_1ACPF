%% ANALYZE_COMMUNITY_NET_POSITION
%
% Per-bus-phase P/Q bookkeeping at the four energy-community nodes (645,
% 611, 652, 671), decomposed and tracked across all three stages
% (S, S', S''), for reuse in later analysis without re-deriving it ad hoc
% each time. Fixes the capacitor-labeling ambiguity found during manual
% verification: 611.c has a local fixed capacitor, so its raw ieee13.m q
% is already net of that capacitor -- this script recovers the TRUE
% customer Q_load explicitly (Q_load_true = q + q_cap) rather than
% reporting the capacitor-netted value as if it were pure load.
%
% Does not touch the power-flow calculations; pure post-processing on
% already-validated pipeline quantities (p, q, p_der, q_der, lambda from
% run_ieee13_pipeline.m).
%
% Outputs:
%   data/community_net_position.mat  (structured, for reuse in Octave/MATLAB)
%   data/community_net_position.csv  (long format, for Excel/Python/etc.)

clear all
close all
clc

run_ieee13_pipeline;
Sbase = 5e6;
phase_names = {'a','b','c'};

rows = struct('bus_id',{},'bus_idx',{},'phase',{}, ...
    'P_load_kW',{},'Q_load_true_kvar',{},'Q_cap_kvar',{},'P_gen_kW',{},'Q_gen_kvar',{}, ...
    'P_net_S_kW',{},'P_net_Sp_kW',{},'P_net_Spp_kW',{}, ...
    'Q_net_S_kvar',{},'Q_net_Sp_kvar',{},'Q_net_Spp_kvar',{});

for k = 1:numel(der_meta.community)
    c = der_meta.community(k);
    for i = 1:numel(c.phase)
        ph = c.phase(i);
        fidx = rw(c.bus, ph);

        r.bus_id  = c.id;
        r.bus_idx = c.bus;
        r.phase   = phase_names{ph};

        r.P_load_kW        = p(fidx)*Sbase/1e3;                    % P has no capacitor component
        r.Q_load_true_kvar = (q(fidx) + q_cap(fidx))*Sbase/1e3;    % recovers TRUE customer Q_load
        r.Q_cap_kvar       = q_cap(fidx)*Sbase/1e3;                 % 0 everywhere except 611.c (100 kvar)
        r.P_gen_kW         = p_der(fidx)*Sbase/1e3;
        r.Q_gen_kvar       = q_der(fidx)*Sbase/1e3;                 % 0 throughout (unity-PF inverters)

        r.P_net_S_kW   = p(fidx)*Sbase/1e3;
        r.P_net_Sp_kW  = p_S1(fidx)*Sbase/1e3;
        r.P_net_Spp_kW = p_S2(fidx)*Sbase/1e3;

        r.Q_net_S_kvar   = q(fidx)*Sbase/1e3;
        r.Q_net_Sp_kvar  = q_S1(fidx)*Sbase/1e3;
        r.Q_net_Spp_kvar = q_S2(fidx)*Sbase/1e3;

        rows(end+1) = r; %#ok<AGROW>
    end
end

% Feeder-wide totals (all 29 non-slack physical terminals), for context.
nonSlack = 4:numel(p);
totals.P_net_S_kW    = sum(p(nonSlack))*Sbase/1e3;
totals.P_net_Sp_kW   = sum(p_S1(nonSlack))*Sbase/1e3;
totals.P_net_Spp_kW  = sum(p_S2(nonSlack))*Sbase/1e3;
totals.Q_net_S_kvar   = sum(q(nonSlack))*Sbase/1e3;
totals.Q_net_Sp_kvar  = sum(q_S1(nonSlack))*Sbase/1e3;
totals.Q_net_Spp_kvar = sum(q_S2(nonSlack))*Sbase/1e3;
totals.lambda = lambda;
totals.total_DER_kW = der_meta.total_DER_kW;

notes = [ ...
    'P_net = P_load - P_gen (net consumption convention, positive = drawn from network). ', ...
    'S: P_gen=0. S'': P_load unchanged from S, P_gen = installed DER. ', ...
    'S'''': P_load = lambda*P_load(S) [DER NOT scaled]; same pattern for Q via Q_load_true/Q_cap. ', ...
    '611.c has a local fixed 100 kvar capacitor -- Q_load_true already backs it out; ', ...
    'Q_net values include the (unscaled-by-lambda) capacitor contribution, matching what the solver sees.'];

printf('\n============================================================\n');
printf('Community net-position table (%d rows) -- see data/community_net_position.csv\n', numel(rows));
printf('============================================================\n');
printf('%-4s %-3s | %8s %8s %8s | %9s %9s %9s\n', ...
    'bus','ph', 'Pnet(S)','Pnet(S'')','Pnet(S'''')', 'Qnet(S)','Qnet(S'')','Qnet(S'''')');
for i = 1:numel(rows)
    r = rows(i);
    printf('%-4d %-3s | %8.1f %8.1f %8.1f | %9.1f %9.1f %9.1f\n', ...
        r.bus_id, r.phase, r.P_net_S_kW, r.P_net_Sp_kW, r.P_net_Spp_kW, ...
        r.Q_net_S_kvar, r.Q_net_Sp_kvar, r.Q_net_Spp_kvar);
end

%% Save

if ~exist('data','dir'); mkdir('data'); end

save(fullfile('data','community_net_position.mat'), 'rows', 'totals', 'notes', '-v7');
fprintf('\nSaved data/community_net_position.mat\n');

csv_path = fullfile('data','community_net_position.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, ['bus_id,phase,P_load_kW,Q_load_true_kvar,Q_cap_kvar,P_gen_kW,Q_gen_kvar,', ...
    'P_net_S_kW,P_net_Sp_kW,P_net_Spp_kW,Q_net_S_kvar,Q_net_Sp_kvar,Q_net_Spp_kvar\n']);
for i = 1:numel(rows)
    r = rows(i);
    fprintf(fid, '%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', ...
        r.bus_id, r.phase, r.P_load_kW, r.Q_load_true_kvar, r.Q_cap_kvar, r.P_gen_kW, r.Q_gen_kvar, ...
        r.P_net_S_kW, r.P_net_Sp_kW, r.P_net_Spp_kW, r.Q_net_S_kvar, r.Q_net_Sp_kvar, r.Q_net_Spp_kvar);
end
fclose(fid);
fprintf('Saved %s\n', csv_path);
