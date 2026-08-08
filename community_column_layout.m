function layout = community_column_layout(community)
%COMMUNITY_COLUMN_LAYOUT Single source of truth for the 14-column ordering
%of the community flexibility vector x_F, matching A_F_sh's ACTUAL columns.
%
%   layout = COMMUNITY_COLUMN_LAYOUT(community)
%
%   Column order is PER-PLAYER, [p_phases, q_phases] contiguous within each
%   player, players concatenated in community order (645, 611, 652, 671):
%
%     [p645b,p645c,q645b,q645c, p611c,q611c, p652a,q652a, p671a,p671b,p671c,q671a,q671b,q671c]
%
%   This MUST be the only place this ordering is derived. Previously,
%   derive_shared_constraint.m (correctly, matching the DSO-pipeline spec's
%   "reorder by participant" requirement), ieee13_local_capability.m, and
%   the AC-validation index mapping in run_feasibility_check.m /
%   run_io_map_accuracy.m each independently re-derived "the" column
%   order -- and three of those four disagreed with the fourth (9 of 14
%   columns mismatched), silently corrupting every downstream application
%   of x_F back onto the network. All four now call this function instead.
%
%   OUTPUT (struct layout, all fields 1 x 14 unless noted)
%   labels     : cell array of names, e.g. 'p645b', 'q671c'
%   bus        : bus index (ieee13() numbering) per column
%   phase      : phase index (1=a,2=b,3=c) per column
%   is_q       : logical, true if the column is a reactive-power variable
%   full_idx   : rw(bus,phase) -- index into any 3n x 1 vector (p or q)
%   player_of_col : which player (1..4, community order) owns this column
%   n          : total columns (14)
%   block_start, block_end : per-player [start,end] column range (4x1 each)

    phase_names = {'a','b','c'};
    nPlayers = numel(community);

    labels = {}; bus = []; phase_ = []; is_q = []; player_of_col = [];
    block_start = zeros(nPlayers,1); block_end = zeros(nPlayers,1);
    cursor = 0;

    for k = 1:nPlayers
        c = community(k);
        php = numel(c.phase);
        block_start(k) = cursor + 1;

        for ph = c.phase
            labels{end+1} = sprintf('p%d%s', c.id, phase_names{ph}); %#ok<AGROW>
            bus(end+1) = c.bus; phase_(end+1) = ph; is_q(end+1) = false; %#ok<AGROW>
            player_of_col(end+1) = k; %#ok<AGROW>
        end
        for ph = c.phase
            labels{end+1} = sprintf('q%d%s', c.id, phase_names{ph}); %#ok<AGROW>
            bus(end+1) = c.bus; phase_(end+1) = ph; is_q(end+1) = true; %#ok<AGROW>
            player_of_col(end+1) = k; %#ok<AGROW>
        end

        cursor = cursor + 2*php;
        block_end(k) = cursor;
    end

    layout.labels = labels;
    layout.bus = bus;
    layout.phase = phase_;
    layout.is_q = logical(is_q);
    layout.full_idx = arrayfun(@(b,p) 3*(b-1)+p, bus, phase_);
    layout.player_of_col = player_of_col;
    layout.n = numel(labels);
    layout.block_start = block_start;
    layout.block_end = block_end;
end
