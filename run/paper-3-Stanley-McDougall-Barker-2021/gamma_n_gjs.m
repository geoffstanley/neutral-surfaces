function [g,dg_lo,dg_hi] = gamma_n_gjs(s,t,p,along,alat)

%%%	GAMMA_N:	Label hydrographic data with neutral density
%%%
%%%	USAGE:		[g,dg_lo,dg_hi] = gamma_n(s,t,p,along,alat)
%%%
%%%	DESCRIPTION:	Label a section of hydrographic data at a specified
%%%			location with neutral density
%%%
%%%	PRECISION:	Double
%%%
%%%	INPUT:		s       matrix of salinity (each column being a cast)
%%%			t       matrix of in-situ temperatures
%%%			p	matrix of pressures
%%%			along	vector of longitudes (0,360)
%%%			alat	vector of latitudes (-90,90)
%%%
%%%			NOTE:	missing values must be denoted by NaN's
%%%
%%%	OUTPUT:		g	matrix of gamma_n values
%%%			dg_lo   matrix of gamma_n lower error estimates
%%%			dg_hi   matrix of gamma_n upper error estimates
%%%
%%%			NOTE:	NaN denotes missing input data
%%%				-99.0 denotes algorithm failed
%%%				-99.1 denotes input data is outside the valid
%%%				      range of the present equation of state
%%%
%%%	UNITS:		salinity	psu (IPSS-78)
%%%			temperature	degrees C (IPTS-68)
%%%			pressure	db
%%%			gamma_n		kg m-3
%%%
%%%
%%%	AUTHOR:		David Jackett
%%%
%%%	CREATED:	October, 1994
%%%
%%%	REVISION:	2.1		16/2/95
%%%
%%%
%%% Modifications by Geoff Stanley.



%%%
%%%		check # arguments and initialize
%%%

if nargin ~= 5
    error('ERROR in gamma_n.m: invalid input arguments')
end

alat = alat(:).';
along = along(:).';

[nz,nx] = size(s);

gooddata = ~isnan(s+t+p);

if nz == 1
    np = gooddata ;
else
    np = sum(gooddata);
end

along = mod(along, 360);

badcolumn = alat < -80 | alat > 64 | np == 0;
s(:, badcolumn)  = [];
t(:, badcolumn)  = [];
p(:, badcolumn)  = [];
along(badcolumn) = [];
alat(badcolumn)  = [];
np(badcolumn) = [];
indgoodcolumn = find(~badcolumn);
nx2 = size(s,2);
gooddata(:, badcolumn) = [];

g = nan(nz,nx);
dg_lo = nan(nz,nx);
dg_hi = nan(nz,nx);


%%%
%%%		save appropriate array
%%%


savearray = zeros(sum(np) + nx2, 3);
irow = 1;
for ix = 1:nx2
    savearray(irow,:) = [along(ix) alat(ix) np(ix)];
    savearray(irow+1:irow+np(ix), 1) = s(gooddata(:,ix),ix);
    savearray(irow+1:irow+np(ix), 2) = t(gooddata(:,ix),ix);
    savearray(irow+1:irow+np(ix), 3) = p(gooddata(:,ix),ix);
    irow = irow + np(ix) + 1;
end


%%% 		find path to Matlab code
gamma_path = fileparts(which('gamma_n'));
zpwd = pwd;
cd(gamma_path);

save([gamma_path, filesep(), 'glabel_matlab.in'], 'savearray', '-ascii');


%%%
%%%		run external code
%%%

eval('!./glabel_matlab');

glabel_matlab = load([gamma_path, filesep(), 'glabel_matlab.out']);

delete([gamma_path, filesep(), 'glabel_matlab.in']);
delete([gamma_path, filesep(), 'glabel_matlab.out']);


cd(zpwd)



%%%
%%%		assemble gamma_n labels
%%%
if ~isempty(glabel_matlab)
    start = 1;
    for ix = 1:nx2
        n = np(ix);
        finish = start+n-1;
        g(gooddata(:,ix),indgoodcolumn(ix))     = glabel_matlab(start:finish,1);
        dg_lo(gooddata(:,ix),indgoodcolumn(ix)) = glabel_matlab(start:finish,2);
        dg_hi(gooddata(:,ix),indgoodcolumn(ix)) = glabel_matlab(start:finish,3);
        start = finish+1;
    end
end