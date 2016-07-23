function [keval,x,area,d] = meanGenKimsCurveX(X,shapes,L,T120)
%%
Ns=length(shapes);
N=length(shapes{1}.X);
fsampl=125;

%% 
cnt=0;
keval=0;
d=0;

for ii=1:(Ns-1)
    eval(sprintf('shape1=shapes{%d};',ii));
    eval(sprintf('Phi1=shapes{%d}.Phi;',ii));
    
    for jj=(ii+1):Ns
        eval(sprintf('shape2=shapes{%d};',jj));
        eval(sprintf('Phi2=shapes{%d}.Phi;',jj));
        %%
        eval(sprintf('L120=L{%d,%d};',ii,jj));
        %%
        eval(sprintf('C=X.X%d * X.X%d'';',ii,jj));
        T12 =  Phi2 * C.' * Phi1' * shape1.A;
        
        %%
        d_=norm(T120-T12,'fro');
        d=d+d_;
        
        %% 
        %% Soft error
        [~,keval_,x] = GeneralKimEval_final( shape1, shape2, shape1.A, shape2.A, T12, L120 );

        keval=keval+keval_;        
        cnt=cnt+1;
    end
end
%%
keval=keval/cnt;
keval=keval/100;
d=d/cnt;
area = areanum( x, keval );
end