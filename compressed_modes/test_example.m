clear
n=512;
nh=7;

[H,tx]=createH(n);
[X0,~] = svd(randn(n,nh),0);

%{
c=cputime;
[XO,XOcost bO] = OSHER(H,nh,50,1,3,200,X0);
e=cputime;
eO=e-c

c=cputime;
[XN,XNcost bN] = NEUMANN(-H,nh,1/50,4,200,X0);
e=cputime;
eN=e-c

c=cputime;
 [XL, XLcost bL]=SL1_Manopt(H,nh,50,10^(-5),100,X0);
e=cputime;
eL=e-c
%}

n=size(H,1);
% N=20;
% [X0,~] = svd(randn(n,N),0);
c=cputime;
[XM,XMcost bM tM] = MADMM_comptr(H,nh,1*(1/50),2,3,3.5e1,X0);
%[XM,XMcost bM tM] = MADMM_comptr(H,nh,1*(1/50),1,3,3.5e1,X0);
%% Trust regions is much faster than Conjugate Gradient for this example
% [XM,XMcost bM tM] = MADMM_compcg(H,nh,1*(1/50),1,3,25.5e1,X0);

e=cputime;
eM=e-c

X=XM;
for i=1:nh
    h=figure;
   plot(tx,X(:,i),'LineWidth',3,'Color','k');
%    title(num2str(i));
   saveas(h,sprintf('cmm1d-madmm-%d.eps',i),'epsc');
   close all;
%    pause;
end