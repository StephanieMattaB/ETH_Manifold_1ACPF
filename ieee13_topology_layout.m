function [xy, edges, bus_labels] = ieee13_topology_layout()
%IEEE13_TOPOLOGY_LAYOUT Hand-placed (x,y) node coordinates and edge list for
%the IEEE-13 feeder, matching the topology decoded from ieee13.m's
%incidence matrix (branch list: 650-632, 632-633, 633-634, 632-645,
%645-646, 632-671, 671-692, 692-675, 671-680, 671-684, 684-611, 684-652).
%
%   [xy, edges, bus_labels] = IEEE13_TOPOLOGY_LAYOUT()
%
%   xy         : 13 x 2, coordinates in bus-index order (matches ieee13()'s
%                bus order: 650,632,633,634,645,646,671,692,675,680,684,611,652)
%   edges      : 12 x 2, [from_idx, to_idx] pairs (bus-index order)
%   bus_labels : 1 x 13 cell, bus names in the same order

    bus_labels = {'650','632','633','634','645','646','671','692','675','680','684','611','652'};

    xy = [ ...
        0.0,  0.0 ; ...  % 650
        1.2,  0.0 ; ...  % 632
        2.4,  1.8 ; ...  % 633
        3.6,  1.8 ; ...  % 634
        2.4, -1.8 ; ...  % 645
        3.6, -1.8 ; ...  % 646
        2.4,  0.0 ; ...  % 671
        3.6,  0.9 ; ...  % 692
        4.8,  0.9 ; ...  % 675
        3.6,  0.0 ; ...  % 680
        3.6, -0.9 ; ...  % 684
        4.8, -1.4 ; ...  % 611
        4.8, -0.3 ];     % 652

    edges = [ ...
        1  2; ...   % 650-632
        2  3; ...   % 632-633
        3  4; ...   % 633-634
        2  5; ...   % 632-645
        5  6; ...   % 645-646
        2  7; ...   % 632-671
        7  8; ...   % 671-692
        8  9; ...   % 692-675
        7 10; ...   % 671-680
        7 11; ...   % 671-684
        11 12; ...  % 684-611
        11 13];     % 684-652
end
