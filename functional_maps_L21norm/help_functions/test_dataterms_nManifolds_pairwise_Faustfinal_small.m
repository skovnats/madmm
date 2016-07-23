%
clear; close all;
%% Script for comparing different the chaing of 3 shapes to the pairs of shapes
k=25; % number of Fourier bases
niter=8e0; % number of inner iterations
initer=4;%4      % number of inner iterations
%%
mu = 5e1; %% L21 norm weight (1/mu)

%%
% Ns=100; % 0-99
inds=0:7:99; % indexes of shapes
inds=0:65:99;
Ns=length(inds);
% Ns=3;
shape=loadfaust_i(0,1);
%% Number of functions
nf=size(shape.F,2);

for ii=1:(Ns)
    i=inds(ii);
    %% number of points of shapei
    eval(sprintf('N%d=length(shape.X);',ii));
    eval(sprintf('shapes{%d}=loadfaust_i(%d,%d);',ii,i,k));
    eval(sprintf('F{%d}=shapes{%d}.F;',ii,ii));
    eval(sprintf('FourCoeffs{%d} = shapes{%d}.fcoef;',ii,ii));
end
%% Define the ground-truth correspondences
L12=[(1:N1)', (1:N1)'];

for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        %% L12
        L{ii,jj} = L12;
    end
end
%


%% Calculate the Spectral decomposition
% opt=giveopt('nastya',0);
% opt.num_evecs=k;
for ii=1:Ns
    %     eval(sprintf('[W,A{%d},Phi{%d},DD]=DLBO(shapes{%d},opt);',ii,ii,ii));
    eval(sprintf('Sigma{%d} = diag(shapes{%d}.DD);',ii,ii));
%     eval(sprintf('shapes{%d}.Phi = shapes{%d}.Phi;',ii,ii));
end


%% initialization
for jj=2:Ns
    eval(sprintf('[X0.X1, X0.X%d] = initialize_fourier( k, F{1}, shapes{1}.Phi, shapes{1}.A, F{%d}, shapes{%d}.Phi, shapes{%d}.A);',jj,jj,jj,jj));
end

%% Optimization
tic; 
% L21
% [X cXY nnr nns rrho times]=procXYmnp21_n_pairwise(FourCoeffs,1/mu,Sigma,X0,niter,initer);
[X cXY nnr nns rrho times]=procXYmnp21_n(FourCoeffs,1/mu,Sigma,X0,niter,initer);
% C12=X.V*X.W';
t=toc;
%%
save('test-res-L21norm-nManifoldsFaustfinal_small','X');
load('test-res-L21norm-nManifoldsFaustfinal_small','X');


%% MAKS's method
for ii=1:(Ns-1)
   for jj=(ii+1):Ns
       eval(sprintf('CL2g{%d,%d}=L2Norm(FourCoeffs{%d},FourCoeffs{%d});',ii,jj,ii,jj));
   end
end
save('test_FaustMaks_small', 'CL2g');
load('test_FaustMaks_small', 'CL2g');



%% Method of Amit
LL=zeros(Ns*k,Ns*k);
I=eye(k,k);
for ii=1:(Ns-1)
   for jj=(ii+1):Ns
%        LL=blkdiag(LL,CL2g{ii,jj});
        LL((k*(ii-1)+1):(k*(ii)),(k*(jj-1)+1):(k*(jj)))=CL2g{ii,jj};
   end
end
LL=LL'+LL;
for ii=1:Ns;
   LL((k*(ii-1)+1):(k*(ii)),(k*(ii-1)+1):(k*(ii)))=I; 
end

% Find first k eigs
LL=0.5*(LL+LL.');
[ evecs, evals ] = giveEigVV( LL, 1, k );
for ii=1:Ns
   VV{ii}=evecs((1+k*(ii-1)):(k*ii),:); 
%    [u,s,v]=svd(VV{ii});
%    VVo{ii}=u*v';
   [VVo{ii},~]=qr(VV{ii});
end
for ii=1:(Ns-1)
   for jj=(ii+1):Ns
       CA{ii,jj}=VV{jj}*VV{ii}';
       CAo{ii,jj}=VVo{jj}*VVo{ii}';
   end
end
save('test_Amit_Faust_small','VV','VVo','CA','CAo');
load('test_Amit_Faust_small','VV','VVo','CA','CAo');


%%
% T = shapes{jj}.Phi * CL2g{ii,jj}' * shapes{ii}.Phi.' shapes{ii}.A;
% [deviation,distribution,x] = GeneralKimEval_final( shape1, shape2, AX, AY, T12, L12gr );

cnt=0;
keval=0;
kevalL2=0;
kevalA=0;
kevalAo=0;

for ii=1:(Ns-1)
    % for ii=1:Ns
    %%
    eval(sprintf('N=N%d;',ii));
    eval(sprintf('shape1=shapes{%d};',ii));
    eval(sprintf('Phi1=shapes{%d}.Phi;',ii));
    for jj=(ii+1):Ns
        eval(sprintf('shape2=shapes{%d};',jj));
        eval(sprintf('Phi2=shapes{%d}.Phi;',jj));
        %%
        eval(sprintf('L120=L{%d,%d};',ii,jj));
        
        %% ours
        eval(sprintf('C=X.X%d * X.X%d'';',ii,jj));
        %       eval(sprintf('A1=X.X%d; A2 = X.X%d'';',ii,jj));
        %% point-wise correspondence
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        %       [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, eye(size(C)), Phi1*A1, Phi2*A2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final2( shape1, shape2, L120, L12(1:25:end,:) );
        keval=keval+keval_;
        
        %% L2
%         eval(sprintf('CL2=L2Norm(FourCoeffs{%d},FourCoeffs{%d});',ii,jj));
        eval(sprintf('CL2=CL2g{%d,%d};',ii,jj));
        %%
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, CL2, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final2( shape1, shape2, L120, L12(1:25:end,:) );
        %%
        kevalL2=kevalL2+keval_;
        %%
        %% Amit
        eval(sprintf('C=CA{%d,%d};',ii,jj));
        %%
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final2( shape1, shape2, L120, L12(1:25:end,:) );
        %%
        kevalA=kevalA+keval_;
        
        eval(sprintf('C=CAo{%d,%d};',ii,jj));
        %%
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final2( shape1, shape2, L120, L12(1:25:end,:) );
        %%
        kevalAo=kevalAo+keval_;
        
        cnt=cnt+1;
    end
end
%%
keval=keval/cnt;
kevalL2=kevalL2/cnt;
kevalA=kevalA/cnt;
kevalAo=kevalAo/cnt;
%%
area = areanum( x, keval );
areaL2 = areanum( x, kevalL2 );
areaA = areanum( x, kevalA );
areaAo = areanum( x, kevalAo );

%% plotting the results
cols=gencols(4);
plot(x,kevalL2,'Color',cols(1,:),'LineWidth',2);hn;
plot(x,keval,'Color',cols(2,:),'LineWidth',2);hn;
plot(x,kevalA,'Color',cols(3,:),'LineWidth',2);hn;
plot(x,kevalAo,'Color',cols(4,:),'LineWidth',2);hn;

legend([' Ovsjanikov et al.'],[' MADMM'], ['Amit'], ['Amito']);
% legend(['L2-' num2str(areaL2)],['MADMM' num2str(area)]);
ylim([0 1]);set(gcf,'Color','w');

save('test_res_Faust_Final_small');

%}

