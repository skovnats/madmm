function [deviation,distribution,x] = GeneralKimEval_final( shape1, shape2, AX, AY, T12, L12gr )
% Generalized Kim's curve evaluation
% d(i\in X) = \sum_{j=1}^{|Y|} d(j,gr(i)) * t_i(j)/sum_s(t_i(s)) * 1/sqrt(A(Y))
% X->1, Y->2, gr(i)-ground-truth corresponding point of point i,
% t_i(j)-function corresponding to delta function d_i (ideally delta function for gr(i)); A(Y)-area of Y
% AX, AY - vectors of local area elements
% According to my estimation, max(deviation)<=diam/AY
% T12=T12-min(T12(:));
% T12=T12/max(T12(:));
% T12=max(T12,0);
T12=abs(T12);
%
thresh=0.65;
dt=1e-4;

%%The deviation will be calculated at predefined verteces of the first
%%shape (fps algorithm)
if isfield(shape1,'idx')
    idx=shape1.idx;
else
    idx=fps(1,150-1,shape1);
end
if isfield(shape2,'pidx')
    pidx=shape2.pidx;
else
    pidx=L12gr(idx,2);
end

if isfield(shape2,'D')
    D=full(shape2.D);
else
    D = fastmarch_idx(shape2,pidx);
end
% if not enough memory
% Consider splitting matrix T into L and R
%
T=T12./repmat( sum(T12), size(T12,1), 1 );%normalizing columns of T, i.e. t_i(y)
AY=sqrt(sum(AY(:)));% AREA
%
T=D*T;%n2 x n1
ind=sub2ind(size(T),1:length(idx),idx(:).');
deviation=T(ind)/AY;
deviation0=zeros(size(shape1.X));
deviation0(idx)=deviation;
%
n1=length(idx);
x=0:dt:thresh;
x(end+1)=x(end)+dt;
stam=hist(deviation,x);
stam=stam(1:(end-1));
x=x(1:(end-1));
distribution=(cumsum(stam)/n1)*100;
end


function [D] = fastmarch_idx(shape,idx)
n=length(shape.X);
D=zeros(length(idx),n);
%
for i = 1:length(idx)
    ii=idx(i);
    %
    d=d_shape2(shape,ii);
    D(i,:)=d(:).';
end
%
d=D(:);
d=d(d~=Inf);
mx=max(d);
D(D==Inf)=mx;
end