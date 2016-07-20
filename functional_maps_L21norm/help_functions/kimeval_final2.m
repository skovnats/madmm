function [ keval, x, dd, rgeod ] = kimeval_final2( shapes, shapet, Lgr, L, opt )
% Kim's evaluation
% assumed that the correspondence from a shape source to a shape target
if ~exist( 'opt', 'var' )
   opt.dstep = 0.005; 
end

%-
if ~isfield( opt, 'dstep' )
    dstep = 0.005;
else
    dstep = opt.dstep; % step for checking how many points are within the given geod error radius
end
%-
if ~isfield( opt, 'rgeod' )
%     sh = shstat( shapet );
%     geoddiam = sh.geoddiam(1);
%     rgeod = sh.geoddiam(1);
%     rgeod = approxdiam(shapet);
%     Area = giveA( shapet );
%%
%     Area= calcArea( shapet );
%     Area=sqrt(sum(Area));
%%
%     Area=sum(Area);
    %
    
    %----------
%     rgeod = 99;
%     Area=200; %% geodesic diameter
    Area = approxgeodradius(shapet);
else
    rgeod = opt.rgeod; % geaod radius of the target shape
end
%-
if ~isfield( opt, 'maxpercentage' )
    maxpercentage = 0.25;
else
    maxpercentage = opt.maxpercentage; % maximal percentage of the error to consider
end
%-

%-
% xyzs = xyzshape( shapes );
% xyzt = xyzshape( shapet );

%- max geod error in number
% maxerror = maxpercentage * rgeod; % 35% of the geodesic error

%- geod errors intervals
% errors_ind = 0:dstep:maxerror; 
errors_ind = 0:dstep:maxpercentage;
errors_ind = errors_ind(:);
errors_ind(end+1)=errors_ind(end)+dstep;
%-
N = size( L, 1 ); % number of points which were matched

%-
tic
dd=zeros(N,1);
for i = 1:N
    inds = L( i, 1 ); %source
    indt = L( i, 2 ); %mapped point
    indtorig = Lgr( inds, 2 );%ground-truth correspondence
    d = d_shape( shapet, indtorig );
    %
    d_dist=d(indt)/Area;%distortion=ditance to the ground-truth
    %
    dd(i)=d_dist;
end

N=length(dd);
% keval = keval/N;
x = errors_ind;
%
stam=hist(dd,x);
% stam=[stam(2),stam(2:(end-1)),stam(end-1)];
keval=(cumsum(stam)/N)*1;
%
x=x(1:end-1);
keval=keval(1:end-1);
toc
end

