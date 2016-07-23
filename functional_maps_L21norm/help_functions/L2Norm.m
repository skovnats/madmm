function [X,Xcost] = L2Norm(A,B)
% MADMM method
% Minimizes s.t. |AX-B|_L2
% on the manifold of n x k- orthogonal matrices.

% INPUT:
% A:        a N x n - matrix
% B:        a N x k - matrix
% k:        number of colums of X


% OUTPUT:
% X is the optimum matrix of the main variable
% Xcost is the optimal value of the cost function
X=A\B;
Xcost=norm(A*X-B,'fro')^2;
end