function out = CC2periodic(CC, per, outtype)
% CC2PERIODIC: Find connected components of periodic binary data.
%   CCper = CC2PERIODIC(CC, per) modifies the connected components CC of a
%   2D binary image BW, determined from CC = bwconncomp(BW), as if BW is
%   periodic in any dimension d for which per(d) is true, when per is a
%   logical vector of length 2, or for which d is an element of per, when
%   per is a vector of length <= 2.
% 
%   L = CC2PERIODIC(CC, per, 'LabelMatrix') returns L, analagous to 
%   labelmatrix(CC), but modified to account for periodicity.
%
% See also BWCONNCOMP, LABELMATRIX
%
% Geoff Stanley
% geoff.stanley@physics.ox.ac.uk
% geoffstanley@gmail.com
% 
% History:
% 05/04/2018 - Enabled support for dimension >= 3 if connectivity is minimal
%            - core algorithm changed, to handle objects that stretch between sides
%            - Updated tests
% 15/02/2018 - documentation update
% 09/05/2017 - initial file
%
% Some examples/tests are given below:
%
% %% Test 2D, 4 connectivity
% conn = 4; % Minimal connectivity for 2D
% per = [0 1]; % Wrap columns
% BW = [ ...
%     1 1 1 0 1; ...
%     0 0 1 0 0; ...
%     1 0 1 1 1; ...
%     0 0 0 0 0; ...
%     1 0 1 1 1; ...
%     0 0 1 0 0; ...
%     1 1 1 0 1; ...
%     0 0 0 0 0; ...
%     1 0 0 1 1; ...
%     0 0 0 1 0; ...
%     1 0 0 1 1; ...
%     0 0 0 0 0; ...
%     1 1 0 0 1; ...
%     0 1 0 0 0; ...
%     1 1 0 0 1;];
% CC = bwconncomp(BW,conn);
% disp('Original labelmatrix:'); 
% disp(labelmatrix(CC))
% disp('After CC2periodic:'); 
% disp(CC2periodic(CC,per,'L'))
% 
% %% Test 2D, 8 connectivity
% conn = 8; % Maximal connectivity for 2D
% per = [0 1]; % Wrap columns
% BW = [ ...
%     1 0 0; ...
%     0 0 1; ...
%     1 0 0; ...
%     0 0 0; ...
%     0 0 1; ...
%     1 0 0; ...
%     0 0 1; ];
% CC = bwconncomp(BW,conn);
% disp('Original labelmatrix:'); 
% disp(labelmatrix(CC))
% disp('After CC2periodic:'); 
% disp(CC2periodic(CC,per,'L'))
% 
% %% Test 3D
% conn = 6; % Minimal connectivity for 3D
% per = [1 0 1]; % Wrap columns and bottom/top
% clear BW
% BW(:,:,1) = [ ...
%     1 0 0; ...
%     0 0 1; ...
%     1 0 0; ];
% BW(:,:,2) = [ ...
%     0 0 0; ...
%     0 1 0; ...
%     0 0 0; ];
% BW(:,:,3) = [ ...
%     1 0 0; ...
%     0 0 1; ...
%     1 1 0; ];
% CC = bwconncomp(BW,conn);
% disp('Original labelmatrix:'); 
% disp(labelmatrix(CC))
% disp('After CC2periodic:'); 
% disp(CC2periodic(CC,per,'L'))

if nargin < 3 || isempty(outtype)
    outtype = 'CC';
end

assert(isstruct(CC) ...
    && isfield(CC, 'Connectivity') ...
    && isfield(CC, 'ImageSize') ...
    && isfield(CC, 'NumObjects') ...
    && isfield(CC, 'PixelIdxList'), ...
    'CC must be a struct, as output from bwconncomp()');

D = length(CC.ImageSize);
assert(D == 2 || CC.Connectivity == D*2, ...
    'Currently only supports 2-D images or N-D images with minimal connectivity');

assert(isvector(per), 'per must be a vector');
if length(per) == D
    per = find(per(:).');
end

if isempty(per) % No periodic dimensions; nothing to do. 
    out = CC;
    return
end

assert(length(per) <= D && min(per) >= 1 && max(per) <= D, ...
    'per must be a vector of length 2 or a vector of length <= D with values 1 and 2.');

assert(ischar(outtype) && upper(outtype(1)) == 'C' || upper(outtype(1)) == 'L', ...
    'outtype must be a string beginning with ''C'' or ''L''');


% Create the matrix giving a unique positive integer to each distinct
% object and having 0 elsewhere.
L = labelmatrix(CC);

if isa(L,'uint8') || isa(L,'uint16') || isa(L,'uint32') || isa(L,'uint64')
    % When L is a uint, negative numbers are automatically converted to 0.
    posdiff = @diff;
else
    % Make diff(L1) >= 0.
    posdiff = @(x) max(0,diff(x));
end

% Create a map from original objects to new objects.
map = 1:CC.NumObjects;

idx = cellstr(repelem(':', D, 1));
for i = per
    idx{i} = 1;
    L1 = L(idx{:}); % e.g. L1 = L(1,:), when D=2 and i=1.
    idx{i} = CC.ImageSize(i);
    L2 = L(idx{:}); % e.g. L2 = L(end,:), when D=2 and i=1.
    idx{i} = ':';
    L1 = L1(:).'; % Ensure row vector, so that for loop with find, below, works
    L2 = L2(:).';
    
    if D == 2 && CC.Connectivity == 8
        % Extend each object backwards one pixel.
        % This enables diagonal connections, but using the simpler code
        % below for non-diagonal connections.
        L1(1:end-1) = L1(1:end-1) + posdiff(L1);
        L2(1:end-1) = L2(1:end-1) + posdiff(L2);
    end
    
    % Find adjacent (non-diagonal) objects on either boundary
    % and pair them up.
    for j = find(L1 ~= L2 & L1 & L2)
        % Collect all original objects that point to the the new object on the 2 side
        % and make them point to the original object on the 1 side.
        map(map == map(L2(j))) = map(L1(j));
    end
    
end

% Collapse labels into 1:length(unique(map))
[mapuniq,~,map] = unique(map);
NumObjects = length(mapuniq);

if upper(outtype(1)) == 'C'
    out = CC;
    out.NumObjects = NumObjects;
    
    % Collect all indices of old objects. 
    idx = accumarray(map, (1:length(map)), [out.NumObjects 1], @(x) {x});

    % Re-form the pixels associated with the new set of objects
    out.PixelIdxList = cell(out.NumObjects,1);
    for i = 1:out.NumObjects
        out.PixelIdxList{i} = sort(vertcat(CC.PixelIdxList{idx{i}}));
    end
else %if upper(outtype(1)) == 'L'
    
    % Determine smallest data type: code lifted from labelmatrix()
    if NumObjects <= intmax('uint8')
        dataType = 'uint8';
    elseif NumObjects <= intmax('uint16')
        dataType = 'uint16';
    elseif NumObjects <= intmax('uint32')
        dataType = 'uint32';
    else
        dataType = 'double';
    end
    out = zeros(CC.ImageSize,dataType);
    
    out(L>0) = map(L(L>0));
end

