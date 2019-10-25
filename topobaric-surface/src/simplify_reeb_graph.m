function                  [arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes] ...
    = simplify_reeb_graph( arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, ...
    nArcRemain, WEIGHT_PERSIST )
%SIMPLIFY_REEB_GRAPH  Simplify the Reeb Graph by leaf pruning.
%
%
% --- Input:
% arc_from, arc_to, arc_id, arc_segment, node_prev, node_next, node_type,
%   node_fn, node_v, nArcs: these define the Reeb graph. For further
%   information, see the documentation in topobaric_surface.
% nArcRemain: [1, 1] number of arcs to remain after simplification (if an
%   integer), or the fraction of the input number of arcs to remain (if a
%   real number between 0 and 1). (There may be a few less than this if
%   more than one arcs are killed on the final step.)
% WEIGHT_PERSIST [1, 1]: relative weight for simplification by persistence,
%   between 0 and 1. The weight of an arc is (WEIGHT_PERSIST) * persistence
%   + (1-WEIGHT_PERSIST) * num, where persistence is the absolute value of
%   the difference between values associated with the two nodes incident to
%   the arc, and num is the total number of vertices associated to the arc.
%
%
% --- Output:
% arc_from, arc_to, arc_id, arc_segment, node_prev, node_next, node_type,
%   node_fn, node_v, nArcs: the simplified Reeb Graph.
%
%
% --- Algorithm:
% This simplification algorithm follows Carr et al. (2010), which is for
% contour trees, i.e. no cycles, but is here applied to graphs with cycles.
% Leaf nodes are iteratively pruned, though the last leaf in either
% direction is protected from pruning. Interior nodes are not changed until
% they become leaf nodes. All cycles in the graph are preserved.
%
% Note: Leaf pruning generally does destroy the property that all vertices
% on an arc have function values in the range of the nodes' function values
% on either end of that arc.
%
% Note: when the domain is not homeomorphic to a disc, there can be (rare)
% nodes with 2 up arcs and no down arcs, or vice versa. (Carr et al. 2010)
%
%
% --- References:
% Carr, H., Snoeyink, J. & van de Panne, M. Flexible isosurfaces:
% Simplifying and displaying scalar topology using the contour tree.
% Computational Geometry 43, 42?58 (2010).

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 2.0.0
%
% Modified by : --
% Date        : --
% Changes     : --

WEIGHT_NUMPIX = 1 - WEIGHT_PERSIST;

if 0 < nArcRemain && nArcRemain < 1
    % Interpret input as a fraction of the original number of arcs
    nArcRemain = round(nArcs * nArcRemain);
end

%nVerts = length(unique(vertcat(arc_segment{:})));

aliveArc = true(nArcs, 1);
nArcAlive = nArcs; % counter for number of arcs remain alive


%% --- Prepare for Leaf Pruning. Build Priority Queue
PQ = PriorityQueue(1);

arc_persist = node_fn(arc_to) - node_fn(arc_from);
WEIGHT_PERSIST = WEIGHT_PERSIST / max(arc_persist);

arc_numpix = cellfun('length', arc_segment);
WEIGHT_NUMPIX = WEIGHT_NUMPIX / max(arc_numpix);

% Weights (per arc) for leaf pruning:
arc_weight = arc_persist * WEIGHT_PERSIST + arc_numpix * WEIGHT_NUMPIX;

% Holes in the manifold (ie. islands in the ocean) were filled in
% pre-processing, then in post-processing any node in such a hole was set
% to have node_v == 0. Set any arc adjacent to such a node to have weight
% 0, so it gets pruned first.
arc_weight(node_v(arc_to) == 0 | node_v(arc_from) == 0) = 0;

degup = cellfun('length', node_next);
degdn = cellfun('length', node_prev); % maxima have node_type == 3
leafmin = find(degup == 1 & degdn == 0);
leafmax = find(degup == 0 & degdn == 1);

% --- Add leaf minima when their parent is dn-forking saddle
for leaf = leafmin(:).'
    arc = node_next{leaf};
    if aliveArc(arc) && degdn(arc_to(arc)) > 1
        PQ.insert([arc_weight(arc), leaf]);
    end
end

% --- Add leaf maxima when their parent is up-forking saddle
for leaf = leafmax(:).'
    arc = node_prev{leaf};
    if aliveArc(arc) && degup(arc_from(arc)) > 1
        PQ.insert([arc_weight(arc), leaf]);
    end
end


% Don't use the following in the pruning loop, because things will have
% changed!
clear degup degdn arc_weight arc_persist arc_numpix


% Arc Partition: arcPart{i} is a vector of all arcs (in RG) that are
% collapsed onto arc i.
arcPart = num2cell((1:nArcs).');

%% Begin leaf pruning:
while (PQ.size() > 0 && nArcAlive > nArcRemain)
    
    % Get lowest persistence:
    pop = PQ.remove();
    w_  = pop(1);
    leaf = pop(2);
    %assert(w_ >= 0, 'Negative priority!');
    
    
    if node_type(leaf) == 1 % minimum
        
        arc = node_next{leaf};
        par = arc_to(arc);          % parent
        oth = node_prev{par};       % other
        
        % Re-check priority. If it's changed, put it back:
        w = getweight(par, leaf, arc);
        if w ~= w_
            PQ.insert([w, leaf]);
            continue
        end
        
        if length(oth) >= 2 && aliveArc(arc)
            % Leaf arc can be merged with arc to other child.
            
            % Prune connection from parent to leaf:
            oth = oth(oth ~= arc);
            node_prev{par} = oth;
            
            if length(oth) >= 2 % Must choose which other to merge with
                
                % Choose the most minimum one
                [~,i] = min(node_fn(arc_from(oth)));
                oth = oth(i);
                
                % Merge arc segments:
                arcPart{oth} = vertcat(arcPart{[oth; arc]});
                arcPart{arc} = [];
                aliveArc(arc) = false;
                nArcAlive = nArcAlive - 1;
                
            else
                % Merge arc segments:
                arcPart{oth} = vertcat(arcPart{[oth; arc]});
                arcPart{arc} = [];
                aliveArc(arc) = false;
                nArcAlive = nArcAlive - 1;
                
                % Guaranteed: length(oth) == 1
                eld = node_next{par};   % elder
                sib = arc_from(oth);    % sibling of leaf node
                sibup = length(node_next{sib}); % Number of arcs up from sibling
                
                if length(eld) == 1
                    % Vertex reduction!
                    
                    % Change gpar <-> par connection to become gpar <-> sib. (delete eld)
                    gpar = arc_to(eld);
                    node_prev{gpar}(node_prev{gpar} == eld) = oth;
                    arc_to(oth) = gpar;
                    %assert(aliveArc(eld), 'This arc already dead!')
                    arcPart{oth} = vertcat(arcPart{[oth; eld]});
                    arcPart{eld} = [];
                    aliveArc(eld) = false;
                    nArcAlive = nArcAlive - 1;
                    
                    % Check if gpar is a prunable leaf
                    prune_up = isempty(node_next{gpar}) && length(node_prev{gpar}) == 1 && sibup > 1;
                    
                    % Check if sibling is prunable leaf
                    prune_dn = isempty(node_prev{sib}) && sibup == 1 && length(node_prev{gpar}) > 1;
                    
                    if prune_up || prune_dn
                        w = getweight(gpar, sib, oth);
                        
                        if prune_up
                            PQ.insert([w, gpar]);
                        end
                        
                        if prune_dn
                            % The priority calculated here is the same as
                            % the priority for the above gpar, as it's the
                            % same arc.
                            PQ.insert([w, sib]);
                        end
                    end
                    
                elseif isempty(eld) && sibup > 1 % Second condition just checks prunable
                    % Rare case: parent had 0 up arcs, 2 down arcs.
                    % par remains the parent of child = arc_from(oth).
                    % Even if leaf < child, leave edge from child to par as an up-arc
                    % (this just includes additional low data in the region of this edge;
                    % doesn't change that to get to leaf from child, you must initially go up).
                    w = getweight(par, sib, oth);
                    PQ.insert([w, par]);
                end
            end
        end
        
    else % should have node_type(leaf) == 3  % maximum
        %assert(node_type(leaf) == 3, 'not a maximum!');
        
        arc = node_prev{leaf};
        par = arc_from(arc);        % parent
        oth = node_next{par};       % other
        
        % Re-check priority. If it's changed, put it back:
        w = getweight(leaf, par, arc);
        if w ~= w_
            PQ.insert([w, leaf]);
            continue
        end
        
        if length(oth) >= 2 && aliveArc(arc)
            % Leaf arc can be merged with arc to other child.
            
            % Prune connection from parent to leaf
            oth = oth(oth ~= arc);
            node_next{par} = oth;
            
            if length(oth) >= 2 % Must choose which other to merge with
                
                % Choose the most maximum one
                [~,i] = max(node_fn(arc_to(oth)));
                oth = oth(i);
                
                % Merge arc segments:
                arcPart{oth} = vertcat(arcPart{[oth; arc]});
                arcPart{arc} = [];
                aliveArc(arc) = false;
                nArcAlive = nArcAlive - 1;
                
            else
                % Merge arc segments:
                arcPart{oth} = vertcat(arcPart{[oth; arc]});
                arcPart{arc} = [];
                aliveArc(arc) = false;
                nArcAlive = nArcAlive - 1;
                
                % Guaranteed: length(oth) == 1
                eld = node_prev{par};   % elder
                sib = arc_to(oth);      % sibling of leaf node
                sibdn = length(node_prev{sib}); % Number of arcs dn from sibling
                
                if length(eld) == 1
                    % Vertex reduction!
                    
                    % Change gpar <-> par connection to become gpar <-> sib. (delete eld)
                    gpar = arc_from(eld);
                    node_next{gpar}(node_next{gpar} == eld) = oth;
                    arc_from(oth) = gpar;
                    %assert(aliveArc(eld), 'This arc already dead!')
                    arcPart{oth} = vertcat(arcPart{[oth; eld]});
                    arcPart{eld} = [];
                    aliveArc(eld) = false;
                    nArcAlive = nArcAlive - 1;
                    
                    % Check if gpar is a prunable leaf
                    prune_dn = isempty(node_prev{gpar}) && length(node_next{gpar}) == 1 && sibdn > 1;
                    
                    % Check if sibling is prunable leaf
                    prune_up = isempty(node_next{sib}) && sibdn == 1 && length(node_next{gpar}) > 1;
                    
                    if prune_dn || prune_up
                        w = getweight(sib, gpar, oth);
                        
                        if prune_dn
                            PQ.insert([w, gpar]);
                        end
                        
                        if prune_up
                            % The priority calculated here is the same as
                            % the priority for the above gpar, as it's the
                            % same arc.
                            PQ.insert([w, sib]);
                        end
                    end
                    
                elseif isempty(eld) && sibdn > 1 % Second condition just checks prunable
                    % Rare case: parent had 0 up arcs, 2 down arcs.
                    % par remains the parent of child = arc_from(oth).
                    % Even if leaf < child, leave edge from child to par as an up-arc
                    % (this just includes additional low data in the region of this edge;
                    % doesn't change that to get to leaf from child, you must initially go up).
                    w = getweight(sib, par, oth);
                    PQ.insert([w, par]);
                end
            end
        end
        
    end
    
    
    %degup = cellfun('length', node_next);
    %degdn = cellfun('length', node_prev);
    %assert( sum(degup==1 & degdn==1) == 0, 'Regular node produced!' );
    
    % Extra slow assertion here:
    %allverts = unique(vertcat(arc_segment{vertcat(arcPart{:})}));
    %assert(length(allverts) == nVerts, 'Lost some vertices!')
    
end

%% Tidy the RG
class_arc = class(arc_from);
arc_from = arc_from(aliveArc);
arc_to   = arc_to(aliveArc);
arc_id   = 1 : nArcs;
arc_id   = arc_id(aliveArc);
nArcs = length(arc_from);

tmp = cell(nArcs,1);
for i = 1:nArcs
    tmp{i} = unique(vertcat(arc_segment{arcPart{arc_id(i)}}));
end
arc_segment = tmp;

% Make a list of the node ID's in the old graph that the RG uses.
[old_nodes,~,iC] = unique([arc_from; arc_to]); % old_nodes(iC) == [arc_from; arc_to]
iC = cast(iC, class_arc); % Convert back to correct class

% Change arcs to go between node indices that are 1 through # of nodes in RG
arc_from = iC(1:end/2);
arc_to = iC(end/2+1 : end);

% These may not accurate for non-saddle nodes! But we don't need them there.
node_v = node_v(old_nodes);
node_type = node_type(old_nodes);
node_fn = node_fn(old_nodes);
nNodes = length(node_v);

% Relabel arcs to be 1 through nArcs, and update node links to point to correct arcs:
node_prev = node_prev(old_nodes);
node_next = node_next(old_nodes);
map_arcid = zeros(length(aliveArc),1,class_arc);
map_arcid(aliveArc) = 1:nArcs;
for n = 1:nNodes
    node_prev{n} = map_arcid(node_prev{n});
    node_next{n} = map_arcid(node_next{n});
    %assert(all(node_prev{n} > zeros(1, 1, 'like', node_prev{1})) && all(node_next{n} > zeros(1, 1, 'like', node_prev{1})), 'Node mapping broken');
end
%assert(~any(isnan(vertcat(node_prev{:}))), 'Node prev links broken!')
%assert(~any(isnan(vertcat(node_next{:}))), 'Node next links broken!')

    function w = getweight(n1, n2, arc)
        
        numpix = sum(cellfun('length', arc_segment(arcPart{arc}))); % Don't worry about occasional duplicate vertex (on nodes)
        if node_v(n1) == 0 || node_v(n2) == 0
            persist = 0;
        else
            persist = node_fn(n1) - node_fn(n2);
        end
        w = persist * WEIGHT_PERSIST + numpix * WEIGHT_NUMPIX;
    end

end
