function [X, Xcost b0 t0]=SL1_Manopt(H,N,mu,eps,it,X0)

% Minimizes |X|_eps/mu+trace(X'HX)
% on the Stiefel manifold X'*X=I

% Here |.|_eps is the smoothed L1- norm |x|=sqrt(x^2+eps), eps>0. 

% INPUT:
% H: discrete Hamiltonian
% N: number of colums of X 
% eps: smoothing parameter for L1 norm (eps = 10^(-6) )
% it: number of iterations of MANOPT conjugate gradient algorithm

% OUTPUT:
% X solution of smoothed L1 minimization by MANOPT
% Xcost optimal value of cost function f(X)=|X|_eps/mu+trace(X'HX)
warning('off');

    n=size(H,1);
    
  % Start matrix  
    %[X0 SX]=eigs(H,N);
        if exist('X0','var')
        X=X0;
    else
        [X,~] = svd(randn(n,N),0);
    end
    %X0=rand(n,N);
    
    b0=cost(H,X0,mu,eps);
    t0=0;

  % Create the problem structure.
    manifold = stiefelfactory(n, N, 1);
    problem.M = manifold;
    
  % Define the problem cost function and its gradient.
    problem.cost = @(X) cost(H,X,mu,eps);
    egrad = @(X) egra(H,X,mu,eps);
    problem.grad = @(Y) manifold.egrad2rgrad(Y, egrad(Y));
    
   % Numerically check the differential
   % checkgradient(problem);

  % Stopfunction
  %  options.stopfun = @mystopfun;
    function stopnow = mystopfun(problem, x, info, last)
      stopnow = (last >= 3 && (info(last-2).cost - info(last).cost)/info(last).cost < 1e-8);
    end
    options.maxiter=it;
%     options.minstepsize=1e-1;
%     options.verbosity=0;
    options.verbosity=2;

  % Solve.
    problem.ff = @(X) trace(X.'*H*X) + sum(abs(X(:)))/mu;
    [X, Xcost, info, options] = trustregions(problem,X0,options);
%     [X, Xcost, info, options] = conjugategradient(problem,X0,options);
  
    Xcost=cost(H,X,mu,eps);   
%     t0=[];
%     b0=[];
    for i=1:size(info,2)
        b0=[b0, info(i).cost];
        t0=[t0, info(i).time];
    end;
end

function cc=cost(H,X,mu,eps);
    % smoothed L1 norm in cost function
    
    cc = sum(sqrt(X(:).^2+eps))/mu+trace(X'*H*X); 

end

function eg=egra(H,X,mu,eps)
% gradient of cost function    
    [n m]=size(X);
      
    g1=X./sqrt(X.*X+eps*ones(n,m));      
    g2=(H+H')*X;

    eg=g1/mu+g2;   
end


