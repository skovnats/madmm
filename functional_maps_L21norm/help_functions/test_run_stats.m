%% Script runs different statistics
load('test_run_methods');

%% Define the ground-truth correspondences | Same for all FAUST shapes
L12=[(1:N1)', (1:N1)'];
for ii=1:(Ns-1)
    for jj=(ii+1):Ns
        %% L12
        L{ii,jj} = L12;
    end
end

%% Hard error - Kim's curves
eval(sprintf('C=X.X%d * X.X%d'';',ii,jj));
[shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
[ keval_, x ] = kimeval_final2( shape1, shape2, L120, L12(1:125:end,:) );
area = areanum( x, keval );
plot(x,kevalL2,'Color',cols(1,:),'LineWidth',2);hn;
ylim([0 1]);set(gcf,'Color','w');

%% Distance to the ground-truth
T12 = shapes{jj}.Phi' * C{ii,jj}.' * shapes{ii}.Phi * shapes{ii}.A;
d=norm(T0-T12,'fro');

%% Soft error
[deviation,distribution,x] = GeneralKimEval_final( shapes{ii}, shapes{jj}, shapes{ii}.A, shapes{jj}.A, T12, L12gr );

%% 