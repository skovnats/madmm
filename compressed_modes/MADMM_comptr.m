function [X,Xcost bm tm] = MADMM_comptr(L,N,lambda,rho,steps,it,X0)

% Manifold ADMM method
% Minimizes lambda*|X|_1+trace(X'LX)
% on the manifold of n x N- orthogonal matrices.

% INPUT:
% L:        is the discretized Hamiltonian, a n x n- matrix
% N:        number of colums of X (approximate eigenvectors)
% lambda:    parameter in cost function
% rho>0:     penalty parameter for ADMM
% steps:      number of inner (Manopt) iterations
% it:       is the number of outer iterations (updating the Lagrange vector U)

% OUTPUT:
% X is the optimum matrix of the main variable
% Xcost is the optimal value of the cost function
% bm: history of function values (outer iteration)

n=size(L,1);

% Initializing all variables
%X=polar_svd(rand(n,N));
%X=rand(n,N);
%[X SX]=eigs(H,N);

if exist('X0','var')
    X=X0;
else
    [X,~] = svd(randn(n,N),0);
end
Z=X;
U=zeros(n,N);

%
bm=lambda*sum(abs(X(:)))+trace(X.'*L*X);
tm=0;

% ADMM outer iteration
t_=cputime;
% t=tic;
for i=1:it
    X=iterX(L,N,X,Z,U,lambda,rho,steps);
    Z=iterZ(X,U,lambda,rho);
    
    U=U+X-Z;
    
    Xcost=lambda*sum(abs(X(:)))+trace(X.'*L*X);
    bm=[bm  Xcost];
%     tm=[tm  toc];
    tm=[tm  cputime-t_];
    fprintf('%d:%f\n',i,Xcost);
end;
% tm=tm-t_;
end

function X=iterX(L,N,X,Z,U,lambda,rho,steps)
n=size(X,1);
%  Create the problem structure.
manifold = stiefelfactory(n, N, 1);
problem.M = manifold;

%  Define the problem cost function and its gradient.
problem.cost = @(X) trace(X'*L*X)+rho*norm(X-Z+U,'fro')^2/2;

egrad = @(X) egra(L,X,Z,U,rho);
problem.grad = @(Y) manifold.egrad2rgrad(Y, egrad(Y));

%  Numerically check the differential
%   checkgradient(problem);

%  Stopfunction
options.stopfun = @mystopfun;
    function stopnow = mystopfun(problem, x, info, last)
        stopnow = (last >= 3 && (info(last-2).cost - info(last).cost)/info(last).cost < 1e-8);
    end
options.maxiter=steps;
options.verbosity=0;

%  Solve.
warning off
[X, Xcost, info, options] = trustregions(problem,X,options);
% [X, Xcost, info, options] = conjugategradient(problem,X,options);
warning on

Xcost=lambda*sum(abs(X(:)))+trace(X'*L*X);

end

function Z=iterZ(X,U,lambda,rho)
%
Z=shrink(X+U,lambda/rho);

    function [z]=shrink(z,l)
        z=sign(z).*max(0,abs(z)-l);
    end
end
function eg=egra(L,X,Z,U,rho)
% gradient of cost function
[n m]=size(X);

g1=(L+L')*X;
g2=rho*(X-Z+U);

eg=g1+g2;
end


