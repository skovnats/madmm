function [opt] = giveopt( varargin )
% function returns the structure opt for calculating the LB

switch nargin
    case 1
        type = varargin{1};
    case 2
        type = varargin{1};
        w = varargin{2};
    case 3
        type = varargin{1};
        w = varargin{2};
        opt = varargin{3};
end

indicator = false;
switch type
    case 'belkin'
        opt.which = 'Art';
        opt.LB_PARAM = 'Belkin';
    case 'gl' % graph Laplacian
        opt.which = 'Art';
        opt.LB_PARAM = 'Belkin';
        indicator = true;
    case 'belkingen'
        opt.which = 'Art';
        opt.LB_PARAM = 'BelkinGen';
    case 'cot'
        opt.which = 'Art';
        opt.LB_PARAM = 'cotW';
    case 'cotbb'
        opt.which = 'BB';
        opt.LB_PARAM = 'cot';
    case 'fem'
        opt.which = 'Dan';
        opt.LB_PARAM = 'FEM';
    case 'nastya'
        opt.which = 'Nastya';
        opt.num_evecs = varargin{2};
        opt.LB_PARAM = 'nastya';
    case 'chung'
        opt.which = 'Chung';
        opt.num_evecs = varargin{2};
    case 'bruno'
        % will be used cotangent scheme proposed by Bruno,
        % which  they derived from Exterior Discrete Calculus
        % paper: Spectral Geometry Processing with Manifold Harmonics
        opt.which = 'Bruno';
        opt.num_evecs = varargin{2};
end

fldnms = fieldnames( opt );

%% opt parameter
% opt - structure of information for Laplacian or LBO.
%{
opt.which =  'Art';% 'Art'/'BB'/'Dan'
if strcmp( opt.which, 'BB' )
    opt.LB_PARAM = 'cot';% 'neu', 'dir', 'cot', 'euc', 'geo' ('BB')
    % 'neu' - 1st order FEM Neumann
    % 'dir' - 1st order FEM Dirichlet
    % 'cot' - cotangent weights
    % 'euc' - euclidean weights
    % 'geo' - geodesic weights
elseif strcmp( opt.which, 'Art' )
    %- last update - 21.12.2010 - 'BelkinGen'
    opt.LB_PARAM = 'BelkinGen'; % 'cotW', 'Belkin', 'GL', 'BelkinGen'  % ('Art')
    %-
    % 'cotW' - cotangent weights, W - additional weight will be assigned to
    % color (all weights will be of form: exp(|color_1 - color_2|)),
    % multiplacitivetely
    % 'Belkin' - will be defined Laplacian - Beltrami using heat kernel, .i.e.
    % exp with soem provided std.
    % 'GL' -  Graph Laplace
    % NOTE: also possible to build some other general Laplace matrix ( with defined connectivities as function on some
    % distance between vertexes)
else
    opt.LB_PARAM = 'FEM'; % ('Dan')
    % In case if was chosen
end
%}

if ~isfield( opt, 'TYPE_AREA_CALC' )
    opt.TYPE_AREA_CALC = 'gen' ; %'gen', 'genMean' ('Art')
end
if ~isfield( opt, 'TYPE_BC' )
    opt.TYPE_BC = 'Dirichlet' ;%'Neumann', 'Dirichlet' ('Art')
end

if ~isfield( opt, 'CON_TYPE' )
    opt.CON_TYPE = '1Ring'; %'NN', 'RADIUS', 'nRing' ('Art')
    % 'NN' - nearest neighbors will be used for defining connectivity
    % 'RADIUS' - neighbors within some radius will be used for defining connectivity
    % Mainly it'll be used for G(raph)L(aplacian) and Belkin /as described
end

if ~isfield( opt, 'CON_TYPE_CALC' )
    opt.CON_TYPE_CALC = 'loop'; % 'all'/'loop' ('Art')
end
if ~isfield( opt, 'CON_USE_GEODESIC' )
    opt.CON_USE_GEODESIC = false; % 0/1 ('Art')
end

if ~isfield( opt, 'CON_RADIUS' )
    opt.CON_RADIUS = 2; % R, radius or number of nearest neighbprs ('Art')
end

if ~isfield( opt, 'CON_RADIUS_SCALE' )
    opt.CON_RADIUS_SCALE = 1.5; % this scalar will be used, to find first
    % scale*Radius nearest neighbors for each vertix  ('Art')
end

if ~isfield( opt, 'CON_SIGMA' )
    opt.CON_SIGMA = 100*.02;% std for exp for connectivity ('Art')
end

if ~isfield( opt, 'CON_METRIC' )
    opt.CON_METRIC = 'exp'; %, 'euc' 'exp', 'eucMltp' , 'eucExp'('Art')
    if strcmp( type, 'cot' )
        opt.CON_METRIC = 'euc';
    end
end

if ~isfield( opt, 'COLOR_W' )
    opt.COLOR_W = 0;%sqrt(0.1); % 4*opt.CON_SIGMA/(1/2);%scalar weight for color ('Art')
    % Pay attention, that were performed changes to code, so Wnew = sqrt(Wold)
    if exist( 'w', 'var' )==1
        opt.COLOR_W = w;
    end
end

%( 'Art' )

if ~isfield( opt, 'METRIC_FOR_COLOR_COORD' )
    % not relevant in case of 'BelkinGen'
    opt.METRIC_FOR_COLOR_COORD = 'metric';% 'metric'/'eqmetric'/'scaleinv'/'affine' ('Belkin'/'FEM') %%
end

if ~isfield( opt, 'METRIC_FOR_GEOM_COORD' )
    % not relevant in case of 'BelkinGen'
    opt.METRIC_FOR_GEOM_COORD = 'metric';% 'metric'/'eqmetric'/'scaleinv'/'affine' ('Belkin'/'FEM') %%
end

% options for smoothing
% opt.IS_SMOOTH_COLOR = true;
opt.IS_SMOOTH_COLOR = false;
% parameters for smoothing the  color coordinates
opt.Nring = 4;
opt.sigma = 8;
opt.type = 'rgb';
%
opt.IS_SMOOTH_METRIC = false;

%- 'BelkinGen' relevant parametrs, which will be ignored  by other methods
% first of all the weight will be taken only with exp

if ~isfield( opt, 'TYPE_OF_GEOM_COORD_DIST_CALC' )
    opt.TYPE_OF_GEOM_COORD_DIST_CALC = 'metric'; % 'euc'/'metric'/'eqmetric'     % 'BelkinGen' ('Art')
end

if ~isfield( opt, 'TYPE_OF_PHOT_COORD_DIST_CALC' )
    opt.TYPE_OF_PHOT_COORD_DIST_CALC = 'eqmetric'; % 'euc'/'metric'/'eqmetric'     % 'BelkinGen' ('Art')
end
% in case if in above methods was taken specified metric involved
% calculation, then the distances will be calculated from metric, that is
% defined as summ of two First Fundamental Forms ( last one corresponds to photometry data, which
% will be scaled with according specified scalar)

if ~isfield( opt, 'TYPE_OF_AREA_CALC' )
    opt.TYPE_OF_AREA_CALC = 'det'; % 'regular'/'det'    % 'BelkinGen' ('Art')
end
% this option specifies how will be calculated the areas of triangles, and
% as consequence - sum of areas of all triangles, that share some vertex
% i.




% this is in some sense adopting parameter t for color coordinates
opt.LOCAL_SCALE = false; % local scaling of color data ('Art')

if ~isfield( opt, 'num_evecs' )
    opt.num_evecs = 200; %number of biggest eigen values to take
end

opt.fem_deg = 1; % degree of FEM method ('BB')

opt.regularize = false; % if to perform regularization process
% it's needed because when adding color often created not wealth
% triangulation('Art')

opt.regul_thresh = 0.12; % values smaller then given thresh
% will be replaced by thresh - received for regular triang. without
% adding color ('Art')

opt.regul_threshW = 10^(-14); % thresh for areas, all value that are smaller then this
% will be replaced by this number.

opt.hard_regularization = false; % will be added to diagonal above thresh for weight matrix.

opt.WhichDescrDist = 'BB'; % 'Art'/'BB' - for descriptors' distances...
% ... will be used or mine or BB code.

opt.WhatKernelToUse = 'mixed'; % 'geom'/'mixed'
% this option tells which kernerl will be used during finding SSBoF.

% Here will be different options, for manipulating color data

opt.IS_PRECLUSTER = false; % if to pre-cluster color, so instead of
% color vector at each point will be taken it's center.

opt.IS_ADOPT_WEIGHT = false; % !only relevant if above TRUE!
% if was chosen pre-clustering, then
% before applying specified above weight for colors /which are centers/
% they will be scaled in order that there mean will be approximately
% equal to mean edge lenght of shapes, and only after that will be applied
% weightening.

if ~isfield( opt, 'fL' )
    % Finctions of color coordinates ( (L,a,b) -> (fL(L),fa(a),fb(b)) )
    opt.fL = @(x) x; % log(x)..
    opt.fa = @(x) x;
    opt.fb = @(x) x;
end
% Strings for color functions, in order to distiguish results of
% different runs
fL = 'x';
fa = 'x';
fb = 'x';
fstring = [ fL, fa, fb ];

if ~isfield( opt, 'WHICH_COLOR_COORD_USED' )
    opt.WHICH_COLOR_COORD_USED = logical([1 1 1]); % logical vector, specifying which color
    % coordinates will be used, for example if was specified [1,0,1] =>
    % [(l,a,b) -> (L,b)] ( 'Art' )
end

if ~isfield( opt, 'WHICH_COLOR_COORD_2WEIGHT' )
    opt.WHICH_COLOR_COORD_2WEIGHT = logical([1 1 1]); % logical vector, which specifies
end
% which color coordinate will be weighted.
% for example - if
%{
opt.WHICH_COLOR_COORD_USED = [1 0 1];
opt.WHICH_COLOR_COORD_2WEIGHT = [0 0 1]; =>
(L,a,b) ->(fL(L),fa(a),fb(b))->(fL(L),fb(b))->(fL(L),Weight*fb(b))
%}
%%
if ~isfield( opt, 'IS_SUBSAMPLE' )
    opt.IS_SUBSAMPLE = false;
end
if ~isfield( opt, 'sihks' )
    opt.sihks = false;
end
%%
if ~isfield( opt, 'TYPE_METRIC_CALC' )
    opt.TYPE_METRIC_CALC = 'new'; % 'new'/'old'    % 'FEM' ('dan')
    % will be used the interpolation of the x,y,z coord., as well the L, a
    % and b data
end

if ~isfield( opt, 'TYPE_METRIC_COMBINING' )
    opt.TYPE_METRIC_COMBINING = 'sum'; % 'sum'/'comb'/'mltpl'    % 'FEM' ('dan')
    % will be used the 'sum', but theoretically may be used also the
    % 'mltlpl' since the multiplication of 2 and more p.s.d. matrices is a
    % p.s.d. matrix
end

%-
if indicator
    opt.LB_PARAM = 'GL';
end