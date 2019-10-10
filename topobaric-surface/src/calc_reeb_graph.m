function [p, ocean, arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes] = calc_Reeb_graph(p, OPTS)
%CALC_REEB_GRAPH  Calculate the Reeb graph
%
%
% [p, ocean, arc_from, arc_to, arc_segment, node_prev, node_next, node_type,
% node_fn, node_v, nArcs, nNodes] = calc_Reeb_graph(p, OPTS)
% performs the following steps:
% 1, fill small holes with high data in the given rectilinear data p (if
%    OPTS.FILL_PIX is positive or OPTS.FILL_IJ is not empty);
% 2, build a simplical mesh given rectilinear data in p;
% 3, select a single connected component of that mesh;
% 4, perturb the data p on that mesh until all values are unique;
% 5, calculate the Reeb graph of p on that mesh, using the ReCon software
%    (Doraiswamy and Natarajan, 2013);
% 6, post-process the Reeb graph to take out any holes that were filled with
%    high data
% 7, simplify the Reeb graph (if it has more arcs than
%    OPTS.SIMPLIFY_ARC_REMAIN).
%
%
% --- Input:
% p [nx, ny]: rectilinear data
% OPTS [struct]: options. See topobaric_surface documentation for details.
%
%
% --- Output:
% p [nx, ny]: same as the input p, but nan where outside the chosen
%             connected component of the simplical mesh
% ocean, arc_from, arc_to, arc_segment, node_prev, node_next, node_type,
%   node_fn, node_v, nArcs, nNodes: fields of the Reeb graph. See
%   topobaric_surface documentation for details.
%
%
% --- Requirements:
% ReCon, including modifications to call pack() and run() - ../recon/
% bwconncomp, labelmatrix - Image Processing Toolbox
% CC2periodic - https://www.mathworks.com/matlabcentral/fileexchange/66079
% latlon2vertsfaces, perturb_until_unique, SimplifyReebGraph
%
%
% --- References:
% Doraiswamy, H. & Natarajan, V. Computing Reeb Graphs as a Union of
% Contour Trees. IEEE Transactions on Visualization and Computer Graphics
% 19, 249?262 (2013).

% --- Copyright:
% Copyright 2019 Geoff Stanley
%
% This file is part of Topobaric Surface.
%
% Topobaric Surface is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% Topobaric Surface is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with Topobaric Surface.  If not, see
% <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 1.0
%
% Modified by : --
% Date        : --
% Changes     : --

% --- Pre-processing 1: Select one connected region of the ocean, and possibly fill in small holes
DO_FILL = OPTS.FILL_PIX > 0 || ~isempty(OPTS.FILL_IJ);
if DO_FILL
    
    Ground = isnan(p);
    
    % Calculate connected components of Ground
    CC_ground = CC2periodic(bwconncomp(Ground, 8), OPTS.WRAP, 'CC');
    
    % Set to false (will not fill) for land regions that have > specified number of pixels:
    if OPTS.FILL_PIX > 0
        
        npix = cellfun('length', CC_ground.PixelIdxList);
        npix = npix(:).';
        for i = find(npix > OPTS.FILL_PIX)
            Ground(CC_ground.PixelIdxList{i}) = false;
        end
        
    end
    
    % Set to false (will not fill) for any land regions that are connected
    % to specified locations (e.g. the continents and major islands)
    if ~isempty(OPTS.FILL_IJ)
        L = labelmatrix(CC_ground);
        
        for I = 1:size(OPTS.FILL_IJ,1)
            i = OPTS.FILL_IJ(I,1);
            j = OPTS.FILL_IJ(I,2);
            if Ground(i,j)
                Ground(CC_ground.PixelIdxList{L(i,j)}) = false ;
            end
        end
    end
    
    % Now fill in: It won't matter what values are put in there, so long as
    % they are all greater than all surrounding pixels.
    fill = find(Ground);
    offset = max(p(:));
    p(fill) = (1+offset) : (length(fill) + offset);
    
end


% --- Pre-processing 2: Make simplical decomposition
% ij2v is a 2D map from pixel space to the vertex list.
% That is, verts(ij2v(i,j),3) = p(i,j).
% Note, nV == sum(ij2v(:) > 0)
[verts,faces,ij2v,nV] = latlon2vertsfaces(p,OPTS.WRAP,[],[],OPTS.DECOMP,OPTS.REF_IJ);
ocean = ij2v > 0; % true where the pixel is in the Simplical Decomposition

% Ensure all vertices have unique values
[verts(:,3), perturbed] = perturb_until_unique(verts(:,3));


% --- Processing: the Reeb Graph (RG)
RAA = vgl.iisc.recon.incore.ReconAlgorithmAug;
faces = int32(faces-1); % -1 because Java wants 0-indexing, but MATLAB does 1-indexing
RAA.run(verts(:,1), verts(:,2), zeros(size(verts,1),1), verts(:,3), ...
    faces(:,1), faces(:,2), faces(:,3));
clear faces

% Bring the RG into MATLAB as a struct.
% First argument (nV): causes vertices numbered above nV to be left out of arc_segment
% These are not only the extra vertices formed from the simplical decomposition, but also
% extra vertices added inside ReCon when splitting the input into loop free regions.
% Second argument (1): shift everything by +1, for MATLAB 1 indexing.
RG = struct(RAA.pack(nV,1));
clear RAA


% Dereference RG for faster execution. Also convert to double (which is not slower overall)
arc_from    = double(RG.arc_from);
arc_to      = double(RG.arc_to);
arc_segment = cellfun(@double, RG.arc_segment, 'UniformOutput', false);
node_fn     = RG.node_fn;
node_next   = cellfun(@double, RG.node_next, 'UniformOutput', false);
node_prev   = cellfun(@double, RG.node_prev, 'UniformOutput', false);
node_type   = RG.node_type; % Leave as int8
node_v      = double(RG.node_v);
nArcs       = RG.nArcs;
nNodes      = RG.nNodes;


% --- Post-processing 1: Update p to match the verts used in the Reeb graph
% (if they were perturbed), and set p to nan on each pixel that is not a
% vertex in the simplical decomposition (i.e. disconnected regions)
if perturbed
    if lower(OPTS.DECOMP(1)) == 'd'
        % The simplical decomposition contains all vertices of the original rectangular grid
        p(ocean) = verts(:,3);
    else
        % Just take those vertices that were from the original rectangular grid
        p(ocean) = verts(1:nV,3);
    end
end
p(~ocean) = nan;
clear verts


% --- Post-processing 2: if holes were filled, fix node_v and arc_segments
if DO_FILL
    % Adjust node_v
    keep_verts = ~Ground(ocean);
    vertmap = cumsum(keep_verts);
    vertmap(~keep_verts) = 0; % Mark, as 0, any vertices (and nodes via node_v) that were inside filled regions
    if lower(OPTS.DECOMP(1)) == 'd'
        node_v = vertmap(node_v);
    else
        node_v(node_v <= nV) = vertmap(node_v(node_v <= nV));
    end
    
    % Remove missing vertices from arc_segment, and update the vertex indices
    for e = 1:nArcs
        arc_segment{e} = vertmap(arc_segment{e}(keep_verts(arc_segment{e})));
    end
    

    % Finally, reset any pixels that were filled to nan
    p(fill) = nan;
    
    % and update ocean
    ocean = isfinite(p);
end

% --- Post-processing 3: Simplify the Reeb Graph:
if OPTS.SIMPLIFY_ARC_REMAIN < nArcs
    assert(OPTS.DECOMP(1) == 'd', 'Reeb graph simplification with CROSS simplical decomposition disabled')
    [arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes] ...
        = SimplifyReebGraph( ...
        arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, ...
        OPTS.SIMPLIFY_ARC_REMAIN, OPTS.SIMPLIFY_WEIGHT_PERSIST ...
        );
end

% --- Integrity checks:
%{
assert(all(ocean(:) == isfinite(p(:))), 'ocean differs from where p is finite');

livenodes = 0 < node_v & node_v <= nV;

p_ = p(ocean);
chk = p_(node_v(livenodes)) == node_fn(livenodes);
assert(all(chk), 'Check failed: node_fn, node_v, and p_ are mismatched.');
clear livenodes chk

for e = 1:nArcs
    assert( sum(e == node_next{arc_from(e)}) == 1, 'Arc from mismatch!' )
    assert( sum(e == node_prev{arc_to(  e)}) == 1, 'Arc to   mismatch!' )
end
for n = 1:nNodes
    for e = node_next{n}.'
        assert( arc_from(e) == n, 'node next broken!');
    end
    for e = node_prev{n}.'
        assert( arc_to(e) == n, 'node prev broken!');
    end
end
%}

