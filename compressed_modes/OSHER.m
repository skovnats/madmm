function [X,Xcost bo, to] = OSHER(H,N,mu,lambda,rho,it,X0);

% Osher's ADMM method
% Minimizes |X|_1/mu+trace(X'HX)
% on the manifold of n x N- matrices.

% INPUT:
% H:        is the discretized Hamiltonian, a n x n- matrix
% N:        number of colums of X (approximate eigenvectors)
% mu:       penalty parameter in cost function
% lambda>0: penalty parameter for ADMM
% rho>0:    second penalty parameter
% steps  NOT USED (is for the inner iteration for minimizing wrt. X)
% it:       is the number of outer iterations (updating the Lagrange vector U)

% OUTPUT:
% X is the optimum matrix of the main variable 
% Xcost is the optimal value of the cost function
% bo: history of function values (outer iteration)

    n=size(H,1);

 % Initializing all variables   
    %X=polar_svd(rand(n,N));
    %X=rand(n,N);
    %[X SX]=eigs(H,N);
    if exist('X0','var')
        X=X0;
    else
        [X,~] = svd(randn(n,N),0);
    end
    P=X;
    Q=X;
    B=zeros(n,N);
    b=zeros(n,N);
    
    %
    bo=sum(abs(X(:)))/mu+trace(X'*H*X);
    to=0;

    
% ADMM outer iteration
% tic;
% to=[];
t_=cputime;
    for i=1:it
        X=iterX(H,P,Q,B,b,lambda,rho);
        P=iterP(X,B);
        Q=iterQ(X,b,mu,lambda);
        
        B=B+X-P;
        b=b+X-Q;
        
        Xcost=sum(abs(X(:)))/mu+trace(X'*H*X);
        bo=[bo  Xcost];
        to=[to, cputime-t_];
%         to=[to, toc];
        %
        fprintf('%d:%f\n',i,Xcost);
    end;
    
end

function X=iterX(H,P,Q,B,b,lambda,rho)
warning('off');

% Dimensions of data
n=size(P,1);

R=rho*(P-B)+lambda*(Q-b);
A=2*H+(lambda+rho)*eye(n,n);
X=linsolve(A,R);
end
 
function Q=iterQ(X,b,mu,lambda)
% This is the shrinking operation on the variable Q
    Q=shrink(X+b,1/(mu*lambda));
    
    function [z]=shrink(z,l)
        z=sign(z).*max(0,abs(z)-l);
    end
end

function P=iterP(X,B)
   [U S V]=svd(X+B,'econ');
   P=U*V';
end

