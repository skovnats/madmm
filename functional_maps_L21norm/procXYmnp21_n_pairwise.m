function [X cXY nnr nns rrho times]=procXYmnp21_n_pairwise(FourCoeffs,mu,Ps,X0,niter,initer)

% Computes the solution of
% min_{Xi'Xi=I} for i=1:n-1 for j=i+1:n
% mu*| FourCoeffs{i}*Xi-FourCoeffs{j}*Xj |_2,1 }
% end end
% + for i=1:n
% + (|off(Xi'*Ps{i}*Xi)|^2
% end
%
%
% Method: MADMM with fixed number of iterations
% Simultaneous minimization on X,Y
%
% INPUT:
% FourCoeffs: are m x n matrices where m>=n (for example: Fourier coefficients)
% mu:  weight factor for first part of cost function. mu>0
% Ps: are n x n- symmetric (for example: diagonal with eigenvalues on
%      diagonal
% X0:  X0.X1 and X0.X2.. are  start matrices n x p where  n >= p. For example,
%      random initialization.
%
% OUTPUT:
% X:    X.Xi are solution matrices, m x n
% cXY:  history of cost function (outer iteration)

%%
if ~exist('niter','var')
    niter=50;
end
if ~exist('initer','var')
    initer=50;
end

[m,n]=size(FourCoeffs{1,1,2});
p=size(X0.X1,2);

% initialization
X=X0;

%% number of shapes
Ns=length(FourCoeffs);

%% initialization of the statistics
for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        eval(sprintf('nnr{%d,%d}=[];',ii,jj));
        eval(sprintf('nns{%d,%d}=[];',ii,jj));
        eval(sprintf('rrho{%d,%d}=[];',ii,jj));
    end
end


%% For each a pair of shapes I need a separate Z variable
for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        eval(sprintf('Z{%d,%d}=FourCoeffs{%d,%d,%d}*X.X%d-FourCoeffs{%d,%d,%d}*X.X%d;',...
            ii,jj,ii,ii,jj,ii,jj,ii,jj,jj));
        eval(sprintf('Zk{%d,%d}=Z{%d,%d};',ii,jj,ii,jj));
        eval(sprintf('U{%d,%d}=ones(m,p);',ii,jj));
        eval(sprintf('rho{%d,%d}=1;',ii,jj)); %initial penalty parameter
    end
end
%% Cost functions
cXY=costfXY(FourCoeffs,mu,Ps,X0,Ns);

%% Define the manifold
for ii=1:(Ns)
    eval(sprintf('Man.X%d=stiefelfactory(n,p);',ii));
end
problem.M = productmanifold(Man);

%% Define the cost


% ITERATION: Solve

times=[0];
tic;
for i=1:niter
    %% opt w.r.t X
    %% X is a structure, which has fields: Xi
    X=iterXY(FourCoeffs,mu,Ps,X,Z,U,rho,problem,initer,Ns,0,0);
    %% upd Z
    for ii=1:(Ns-1)
        for jj=(ii+1):Ns
            eval(sprintf('Z{%d,%d}=iterZ(FourCoeffs{%d,%d,%d},FourCoeffs{%d,%d,%d},X.X%d,X.X%d,U{%d,%d},rho{%d,%d},mu,p);',...
                ii,jj,ii,ii,jj,jj,ii,jj,ii,jj,ii,jj,ii,jj));
        end
    end
    %% upd U
    for ii=1:(Ns-1)
        for jj=(ii+1):Ns
            eval(sprintf('U{%d,%d}=U{%d,%d}+FourCoeffs{%d,%d,%d}*X.X%d-FourCoeffs{%d,%d,%d}*X.X%d-Z{%d,%d};',...
                ii,jj,ii,jj,ii,ii,jj,ii,jj,ii,jj,jj,ii,jj));
        end
    end
    %% updating the rho penalty
    for ii=1:(Ns-1)
        for jj=(ii+1):Ns
            eval(sprintf('R=FourCoeffs{%d,%d,%d}*X.X%d-FourCoeffs{%d,%d,%d}*X.X%d-Z{%d,%d};',ii,ii,jj,ii,jj,ii,jj,jj,ii,jj));
            eval(sprintf('S=rho{%d,%d}*([FourCoeffs{%d,%d,%d}'';-FourCoeffs{%d,%d,%d}'']*(Z{%d,%d}-Zk{%d,%d}));',ii,jj,ii,ii,jj,jj,ii,jj,ii,jj,ii,jj));
            %%
              nr=norm(R,'fro');
              ns=norm(S,'fro');
              %%
              if nr>=10*ns
                  eval(sprintf('rho{%d,%d}=2*rho{%d,%d};',ii,jj,ii,jj));
                  eval(sprintf('U{%d,%d}=U{%d,%d}/2;',ii,jj,ii,jj));
              end
              if ns > 10*nr
                  eval(sprintf('rho{%d,%d}=rho{%d,%d}/2;',ii,jj,ii,jj));
                  eval(sprintf('U{%d,%d}=2*U{%d,%d};',ii,jj,ii,jj));
              end
              %% saving the stats
              eval(sprintf('nnr{%d,%d}=[nnr{%d,%d};nr];',ii,jj,ii,jj));
              eval(sprintf('nns{%d,%d}=[nns{%d,%d};ns];',ii,jj,ii,jj));
              eval(sprintf('rrho{%d,%d}=[rrho{%d,%d};rho{%d,%d}];',ii,jj,ii,jj,ii,jj));    
        end
    end
    %% upd variables Zij
    for ii=1:(Ns-1)
        for jj=(ii+1):Ns
               % The "previous" Z matrix
               eval(sprintf('Zk{%d,%d}=Z{%d,%d};',ii,jj,ii,jj));
        end
    end
    
    %% SAVE THE COST
    cXY=[cXY;costfXY(FourCoeffs,mu,Ps,X,Ns)];
    
    %%
    fprintf('%d:%f\n',i,cXY(end));
    %%
    times=[times;toc];
end
end

function X=iterXY(FourCoeffs,mu,Ps,X,Z,U,rho,problem,initer,ns,gV,gW)
% minimizes  sum_ij mu| Zij |_2,1 + rho/2*| FourCoeffiXi-FourCoeffjXj-Zij+Uij |^2+
%   + sum_i norm(off(Xi'Ps{i}Xi),'fro')^2
% with respect to Xi'*Xi=I
% by performing some steps with Manopt minimization

%% Define the options
% Stopfunction
options.stopfun = @mystopfun;
    function stopnow = mystopfun(problem, x, info, last)
        stopnow = (last >= 3 && (info(last-2).cost - info(last).cost)/info(last).cost < 1e-3);
    end
options.maxiter=initer; % 50 -> Klaus parameter

% Switch on/off details of inner iteration
options.verbosity=0;

%  Define the problem cost function and its gradient.
problem.cost  =  @cost ;

    function f=cost(X)
        f=0;
        for ii=1:(ns-1)
            for jj=(ii+1):ns
                eval(sprintf('f=f+0.5*rho{%d,%d}*norm(FourCoeffs{%d,%d,%d}*X.X%d-FourCoeffs{%d,%d,%d}*X.X%d-Z{%d,%d}+U{%d,%d},''fro'')^2;',...
                    ii,jj,ii,ii,jj,ii,jj,ii,jj,jj,ii,jj,ii,jj));
            end
        end
        for ii=1:(ns)
            eval(sprintf('f=f+norm(off(X.X%d''*Ps{%d}*X.X%d),''fro'')^2;',ii,ii,ii));
        end
    end

problem.grad = @(X) problem.M.egrad2rgrad(X, egrad_n(X));
function g=egrad_n(X)
%         F=B*X.W+Z-U;
%         G=A*X.V-Z+U;
%         
%         g.V=zeros(n,p);
%         g.V=g.V+rho*A'*(A*X.V-F);
%         g.V=g.V+4*P*X.V*off(X.V'*P*X.V);
%         
%         g.W=zeros(n,p);
%         g.W=g.W+rho*B'*(B*X.W-G);
%         g.W=g.W+4*Q*X.W*off(X.W'*Q*X.W);
    
   %% Gradient of the off-diag part
   for ii=1:(ns)
       eval(sprintf('g.X%d=4*Ps{%d}*X.X%d*off(X.X%d''*Ps{%d}*X.X%d);',ii,ii,ii,ii,ii,ii));
   end
   %% Graident on the rho/2|AXi-BXj-Zij+Uij| part
   for ii=1:(ns-1)
        for jj=(ii+1):ns
            eval(sprintf('[gV,gW]=egrad(FourCoeffs{%d,%d,%d},Ps{%d},X.X%d,FourCoeffs{%d,%d,%d},Ps{%d},X.X%d,Z{%d,%d},U{%d,%d},rho{%d,%d});',...
                ii,ii,jj,ii,ii,jj,ii,jj,jj,jj,ii,jj,ii,jj,ii,jj));
            eval(sprintf('g.X%d=g.X%d+gV;',ii,ii));
            eval(sprintf('g.X%d=g.X%d+gW;',jj,jj));
        end
   end
%     
end

%  Numerically check the differential
%  checkgradient(problem); 
%  pause;


%}
% Solve.
%[X, Xcost, info, options] = trustregions(problem,X0,options);
[X, Xcost, infoS, options] = conjugategradient(problem,X,options);

end

function Z=iterZ(A,B,V,W,U,rho,mu,p)

% minimizes  mu*|| Z ||_2,1 + rho/2*|| -A*X.V+B*X.W+Z-U ||^2
% with respect to Z by "shrinking"

F=A*V-B*W+U;

for j=1:p
    Z(:,j)=Shrink(F(:,j),mu/rho);
end

    function [z]=Shrink(z,l)
        nz=norm(z);
        z=max(nz-l,0)*z/nz;
    end
end

function Y=off(X)
Y=X-diag(diag(X));
end

function f=costfXY(FourCoeffs,mu,Ps,X,ns)
f=0;
for ii=1:(ns-1)
    eval(sprintf('V=X.X%d;',ii));
    for jj=(ii+1):ns
        A=FourCoeffs{ii,ii,jj};
        B=FourCoeffs{jj,ii,jj};
        eval(sprintf('W=X.X%d;',jj));
        C=A*V-B*W;
        c=sqrt(sum(C.*C));
        f=f+mu*sum(c);
    end
end
for ii=1:(ns)
%     eval(sprintf('V=X.X%d;',ii));
    eval(sprintf('f=f+norm(off(X.X%d''*Ps{%d}*X.X%d),''fro'')^2;',ii,ii,ii));
end
end

function [gV,gW]=egrad(A,P,V,B,Q,W,Z,U,rho)
F=B*W+Z-U;
G=A*V-Z+U;

% gV=zeros(n,p);
gV=rho*A'*(A*V-F);
% gV=gV+4*P*V*off(V'*P*V);

% gW=zeros(n,p);
gW=rho*B'*(B*W-G);
% gW=gW+4*Q*W*off(W'*Q*W);
end