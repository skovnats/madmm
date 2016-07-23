%
clear; close all;
%% Script for calculating correspondence between n shapes
%% you have to provide correspondences between all pairs of these shapes: (2 choose from n) options
nf=25; % number of functions
k=25; % number of Fourier bases

%% MADMM - params
niter=8e2;       % number of inner iterations
initer=4;%4      % number of inner iterations
mu = 5e1; %% L21 norm weight (1/mu)

%%
noisepower=50; % the geodesic radius error
npout=4; % outliers: this number of points will be altered
outpairs={[1 3],[1 4],[1 5]}; % these pairs will be altered (the outliers will be added to the second shape)

%% load data
load('data');

Ns=length(shapes); %%
for ii=1:(Ns)
    %%
    eval(sprintf('N%d=length(shapes{%d}.X);',ii,ii));
end

%% Define the ground-truth correspondences
for ii=1:(Ns-1)
    for jj=(ii+1):(Ns)
        if jj == Ns
            L{ii,jj}=L12;
        else
            L{ii,jj}=[L120(:,1),L120(:,1)];
        end
    end
end

%% Sampling points
ind=fps(1,nf-1,shapes{1}); ind=ind(:);

%% the corresponding points
ind1=ind(:);
for ii=2:Ns
    eval(sprintf('ind%d=L{1,%d}(ind1,2);',ii,ii));
    eval(sprintf('ind%d0=L{1,%d}(ind1,2);',ii,ii));
end

%% Defining pairwise functions
for ii=1:(Ns-1)
    eval(sprintf('indii=ind%d;',ii));
    for jj=(ii):Ns
        eval(sprintf('indjj=ind%d;',jj));
        %%
        eval(sprintf('F{%d,%d,%d}=full(deltai(N%d, indii));',ii,ii,jj,ii));
        eval(sprintf('F{%d,%d,%d}=full(deltai(N%d, indjj));',jj,ii,jj,jj));
    end
end

%% Altering the correspondences, i.e. introducing outliers
for cnt=1:length(outpairs)
    %%
    ii=outpairs{cnt}(1);
    jj=outpairs{cnt}(2);
    %% original indeces
    eval(sprintf('indii=ind%d;',ii));
    eval(sprintf('indjj=ind%d;',jj));
    %% now need to randomly pick npout points
    if npout
        indout_=randperm(nf);
        indout_=indout_(1:npout);
        %% index of points
        indp=indjj(indout_);
        %%
        [indpresampl] = rsamples( shapes{jj}, indp, noisepower );
        %%
        indjj(indout_)=indpresampl;
    end
    
    %% Changing the corresponding functions
    eval(sprintf('F{%d,%d,%d}=full(deltai(N%d, indjj));',jj,ii,jj,jj));
end


%% Calculate the Spectral decomposition
%% Here you provide calculation of Laplacian
opt=giveopt('nastya',0);
opt.num_evecs=k;
for ii=1:Ns
    eval(sprintf('[W,A{%d},Phi{%d},DD]=DLBO(shapes{%d},opt);',ii,ii,ii));
    eval(sprintf('Sigma{%d} = diag(DD(1:%d));',ii,k));
end

for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        %% Generate the Fouries coefficients
        eval(sprintf('FourCoeffs{%d,%d,%d} = F{%d,%d,%d}.'' * A{%d} * Phi{%d};',ii,ii,jj,ii,ii,jj,ii,ii));
        eval(sprintf('FourCoeffs{%d,%d,%d} = F{%d,%d,%d}.'' * A{%d} * Phi{%d};',jj,ii,jj,jj,ii,jj,jj,jj));
    end
end

%% initialization
for jj=2:Ns
    eval(sprintf('[X0.X1, X0.X%d] = initialize_fourier( k, F{1,1,%d}, Phi{1}, A{1}, F{%d,1,%d}, Phi{%d}, A{%d});',jj,jj,jj,jj,jj,jj));
end

%% Optimization
[X cXY nnr nns rrho times]=procXYmnp21_n_pairwise(FourCoeffs,1/mu,Sigma,X0,niter,initer);
% C12=X.V*X.W';

%% Now for each pair of shapes find the correspondence
cnt=0;
keval=0;
kevalL2=0;
for ii=1:(Ns-1)
    %%
    eval(sprintf('N=N%d;',ii));
    eval(sprintf('shape1=shapes{%d};',ii));
    eval(sprintf('Phi1=Phi{%d};',ii));
    for jj=(ii+1):Ns
        %%
        eval(sprintf('shape2=shapes{%d};',jj));
        eval(sprintf('Phi2=Phi{%d};',jj));
        %%
        eval(sprintf('L120=L{%d,%d};',ii,jj));
        
        %% ours
        eval(sprintf('C=X.X%d * X.X%d'';',ii,jj));
        %       eval(sprintf('A1=X.X%d; A2 = X.X%d'';',ii,jj));
        %% point-wise correspondence
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        %       [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, eye(size(C)), Phi1*A1, Phi2*A2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final( shape1, shape2, L120, L12(1:25:end,:) );
        keval=keval+keval_;
        
        %% L2
        eval(sprintf('CL2=L2Norm(FourCoeffs{%d,%d,%d},FourCoeffs{%d,%d,%d});',ii,ii,jj,jj,ii,jj));
        %%
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, CL2, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final( shape1, shape2, L120, L12(1:25:end,:) );
        %%
        kevalL2=kevalL2+keval_;
        %%
        cnt=cnt+1;
    end
end
%%
keval=keval/cnt;
kevalL2=kevalL2/cnt;
%%
area = areanum( x, keval );
areaL2 = areanum( x, kevalL2 );

%% plotting the results
cols=gencols(2);
plot(x,kevalL2,'Color',cols(1,:),'LineWidth',2);hold on;
plot(x,keval,'Color',cols(2,:),'LineWidth',2);hold on;

legend(['L2'],[' MADMM']);
ylim([0 1]);set(gcf,'Color','w');
%}