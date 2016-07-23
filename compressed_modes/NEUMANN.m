function [X,Xcost bo to] = NEUMANN(L,N,lambda,rho,it,X0);

% Neumann's ADMM method
% Minimizes lambda*|X|_1+trace(X'LX)
% on the manifold of n x N- orthogonal matrices.

% INPUT:
% L:        is the discretized Hamiltonian, a n x n- matrix
% N:        number of colums of X (approximate eigenvectors)
% lambda:    parameter in cost function
% rho>0:     penalty parameter for ADMM
% steps  NOT USED (it is for the inner iteration for minimizing wrt. E)
%        IN OUR PROGRAM WE SOLVE THE LIN. EQU. FOR E EXACTLY (USING
%        LINSOLVE()).
% it:       is the number of outer iterations (updating the Lagrange vector U)

% OUTPUT:
% X is the optimum matrix of the main variable
% Xcost is the optimal value of the cost function
% bo: history of function values (outer iteration)

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
E=X;
S=X;
Ue=zeros(n,N);
Us=zeros(n,N);

%
bo=lambda*sum(abs(X0(:)))-trace(X0'*L*X0);
to=0;

% ADMM outer iteration
t_=cputime;
% tic;
for i=1:it
    X=iterX(L,S,E,Us,Ue);
    E=iterE(L,X,Ue,rho);
    S=iterS(X,Us,lambda,rho);
    
    Ue=Ue+X-E;
    Us=Us+X-S;
    
    Xcost=lambda*sum(abs(X(:)))-trace(X'*L*X);
    bo=[bo  Xcost];
    to=[to  cputime-t_];
%     to=[to  toc];
    fprintf('%d:%f\n',i,Xcost);
end;

end

function X=iterX(L,S,E,Us,Ue,rho) % formula (14)

n=size(E,1);
Y=(S-Us+E-Ue)/2;      %
% try
    [U S V]=svd(Y,'econ');
% catch
%     save('Y','Y');
%     pause;
% end
X=U*V';
end

function E=iterE(L,X,Ue,rho)  % formula (17)
n=size(L,1);

R=rho*(X+Ue);
A=(rho*eye(n)-L-L');

E=linsolve(A,R);
% try
    [U S V]=svd(E,'econ');
% catch
%     save('E','E');
%     pause;
% end
end

function S=iterS(X,Us,lambda,rho)
%
S=shrink(X+Us,lambda/rho);

    function [z]=shrink(z,l)
        z=sign(z).*max(0,abs(z)-l);
    end
end