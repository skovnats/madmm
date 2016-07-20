function [varargout] = dstcalc(method, varargin)
% functioon does distance calculations
% inputs: 
% s = dstcalc( 'init', 'diffusion', shape, t, opt );
% d = dstcalc( 'compute', s, sample, shape );

switch lower(method)
    
    case 'init'
        type = varargin{1};
        s = init_distance(type, varargin(2:end));
        varargout(1) = { s };
        
    case 'compute'
        s = varargin{1};
        d = compute_distance(s, varargin(2:end));
        varargout(1) = { d };
        
    case 'deinit'
        s = varargin{1};
        deinit_distance(s);
        
end

function [evals evecs areas W] = compute_spectrum(surface, Nvec)

% Compute Laplacian matrices
opt = giveopt( 'nastya', Nvec );
opt.COLOR_W = 0;
[ W, Am, evecs, evals ] = DLBO( surface, opt );
areas = full( diag( Am ) );

nrm = sqrt(areas(:)'*evecs.^2);
evecs = bsxfun(@rdivide, evecs, nrm); % normalizing evecss


function s = init_distance(type, varargin)

if isfield(varargin{1}{1}, 'type'),
    s = varargin{1}{1};
else
    s = [];
    % here need to add an option of adding the eigendata
    if length( varargin{1} ) > 2
        s  = varargin{1}{3};
    end
    surface = varargin{1}{1};
    s.Nv = length(surface.X);
end
s.type = type;

switch lower(type)
    
    case 'geodesic'
%         s.handle = fastmarchmex('init', int32(surface.TRIV-1), double(surface.X(:)), double(surface.Y(:)), double(surface.Z(:)));
%         s.handle = fastmarchmex('init', int32(surface.TRIV-1), double(surface.X(:)), double(surface.Y(:)), double(surface.Z(:)));
    case 'diffusion'
        t = varargin{1}{2};
        s.kernel = @(lambda)(exp(-t*lambda));
        if ~isfield(s, 'evecs')
            [s.evals, s.evecs s.areas] = compute_spectrum(surface, 200);
        end
        
    case 'commute'
        s.kernel = @(lambda)(1./sqrt(lambda));
        if ~isfield(s, 'evecs')
            [s.evals, s.evecs s.areas] = compute_spectrum(surface, 200);
        end
        
end


function d = compute_distance(s, varargin)

sample = varargin{1}{1};

switch lower(s.type)
    
    case 'geodesic'
        shape=varargin{1}{2};
        %{
        source = repmat(Inf, [s.Nv 1]);
        source(sample) = 0;
        
        
        % TODO: this is extremely inefficient as the grid is re-initialized at
        % each iteration!
        d = fastmarch(shape.TRIV, shape.X, shape.Y, shape.Z, double(source), struct('mode', 'single'));
        
        %         d = fastmarchmex('march', s.handle, double(source));
        d(d>=9999999) = Inf;
        %}
        %-
        d = d_shape( shape, sample );
        
    case {'diffusion', 'commute'}
        kernel = s.kernel;
        d = (bsxfun(@minus, s.evecs(sample,2:end), s.evecs(:,2:end)).^2)*(kernel(s.evals(2:end)).^2);
        d = sqrt(d);
        
    case  'euc'
        shape=varargin{1}{2};
        
        % TODO: this is extremely inefficient as the grid is re-initialized at
        % each iteration!
        d = sqrt((shape.X( sample ) - shape.X).^2+(shape.Y( sample ) - shape.Y).^2+(shape.Z( sample ) - shape.Z).^2);
        
        %         d = fastmarchmex('march', s.handle, double(source));
        d(d>=9999999) = Inf;
end

function deinit_distance(s)
switch lower(s.type)
    
    case 'geodesic'
        fastmarchmex('deinit', s.handle);
        
    case {'diffusion', 'commute'}
        
end



