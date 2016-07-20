function [ W, A, evecs, evals ] = DLBO( shape, opt )
% D(iscrete)L(aplace) B(eltrami) O(perator)
% L = (1./A)*W = inv(A)*W
% shape - structure of shape, it has such 'information' - X, Y, Z, TRIV, in
% case if in metric will be used for using color, then it'll have
% information of color, stored in such manner: shape.R, shape.G...,
% also it can be stored in manner, of Lab color representation(then it'll
% be also stored as shape.R, shape.G..).

% The LBO will be implemented in such maner, that it'll be not important
% wich dimension points are (3-D or 6-D)

% opt - structure of information for Laplacian or LBO.
% opt.which - 'Art'/'BB'/'Dan'/'Nastya'/'Chung'/'Bruno'
% opt.LB_PARAM - 'euc', 'neu', 'dir', 'cot', 'euc', 'geo' ('BB')

% 'neu' - 1st order FEM Neumann
% 'dir' - 1st order FEM Dirichlet
% 'cot' - cotangent weights
% 'euc' - euclidean weights
% 'geo' - geodesic weights

% opt.LB_PARAM - 'cotW', 'Belkin', 'GL' ('Art') + 'FEM' ('Dan')

%  opt.METRIC_FOR_COLOR_COORD  'euc'/'metric'/'eqaff-metric' ('Belkin'/'FEM')

% 'cotW' - cotangent weights, W - additional weight will be assigned to
% color (all weights will be of form: exp(|color_1 - color_2|)),
% multiplacitivetely
% 'Belkin' - will be defined Laplacian - Beltrami using heat kernel, .i.e.
% exp with soem provided std.
% NOTE: also possible to build some other general Laplace matrix ( with defined connectivities as function on some
% distance between vertexes)


% opt.TYPE_AREA_CALC - 'gen',  'genMean' ('Art')

% opt.TYPE_BC - 'Neumann', 'Dirichlet' ('Art')

% opt.CON_TYPE - 'NN', 'RADIUS' ('Art')

% 'NN' - nearest neighbors will be used for defining connectivity
% 'RADIUS' - neighbors within some radius will be used for defining connectivity
% Mainly it'll be used for G(raph)L(aplacian) and Belkin /as described by Mahmudi et al./

% opt.CON_TYPE_CALC - 'all'/'loop' ('Art')

% opt.CON_USE_GEODESIC - 0/1 ('Art')

% opt.CON_RADIUS - R, radius or number of nearest neighbprs ('Art')

% opt.CON_RADIUS_SCALE - this scalar will be used, to find first
% scale*Radius nearest neighbors for each vertix  ('Art')

% opt.CON_SIGMA - sigma, std for exp for connectivity ('Art')

% opt.CON_METRIC - 'euc', 'exp', 'eucMltp' ('Art')

% opt.COLOR_W - scalar weight for color ('Art')

% opt.num_evecs - number of biggest eigen values to take

% opt.fem_deg - degree of FEM method ('BB')

% UPD(31.05.2011) - added additional schemes (Nastya, Chung, Bruno etc)

switch opt.which
    case 'BB'
        switch(opt.LB_PARAM)
            case('cot')
                [evecs, evals, W, A] = main_mshlp('cotangent', shape, opt.num_evecs);
            case('euc')
                [evecs, evals, W, A] = main_mshlp('euclidean', shape, opt.num_evecs);
            case('geo')
                [evecs, evals, W, A] = main_mshlp('geodesic', shape, opt.num_evecs);
            case('neu')
                [evecs, evals, W, A] = fem(shape, 'neu', opt.num_evecs, opt.fem_deg);
            case('dir')
                [evecs, evals, W, A] = fem(shape, 'dir', opt.num_evecs, opt.fem_deg);
            otherwise
                [evecs, evals, W, A] = main_mshlp('cotangent', shape, opt.num_evecs);
                %                 assert(0);
        end
    case 'Art'
        W = giveW( shape, opt );
        if ~strcmp( opt.LB_PARAM, 'GL' )
            A = giveA( shape, opt );
        else
            A = sum( W, 2 );
            A = sparse( 1:length(A), 1:length(A), A );
            %-
            if opt.hard_regularization
                A = full(diag( A )) + opt.regul_threshW;
                A = sparse( 1:length(A), 1:length(A), A );
            end
            W = A - W;
        end
        %          [ evecs, evals ] = giveEigVV( W, A, opt.num_evecs );
    case 'Dan'
        % here photometry can be added only via metric
        [W A] = build_equi_fem_art(shape, opt);
        W = -W;
        %         [ evecs, evals ] = giveEigVV( W, A, opt.num_evecs );
    case 'Nastya'
        [W, A] = calcLB(shape);
        if opt.hard_regularization
            A = max(full(diag( A )),opt.regul_threshW);
            A = sparse( 1:length(A), 1:length(A), A );
        end
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         A = speye( length(A) );
        %!!!!!!!!!!!!!!!!!!!%!!!!!!!!!!!!!!!!!!!
                        W = -W;
        %         [ evecs, evals ] = giveEigVV( W, A, opt.num_evecs );
    case 'Chung'
        [A, W] = FEMLB(shape);
        W = -W;
        %         W=-W;
        %         [ evecs, evals ] = giveEigVV( W, A, opt.num_evecs );
    case 'Bruno'
        [W, A] = calcLBBruno(shape);
        % W = -W; %
        A = 1;
        %         [ evecs, evals ] = giveEigVV( W, 1, opt.num_evecs );
end

if ~strcmp( opt.LB_PARAM, 'GL' )
    W = -W;
end
if nargout > 2
    %     [ evecs, evals ] = giveEigVV( W, 1, opt.num_evecs );
    %        A = 1;
    [ evecs, evals ] = giveEigVV( W, A, opt.num_evecs );
else
    evecs =0 ;
    evals = 0;
end