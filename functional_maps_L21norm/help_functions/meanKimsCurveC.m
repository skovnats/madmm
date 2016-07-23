function [keval,x,area] = meanKimsCurveC(CC,shapes,L)
%%
Ns=length(shapes);
N=length(shapes{1}.X);
fsampl=125;

%% 
cnt=0;
keval=0;

for ii=1:(Ns-1)
    eval(sprintf('shape1=shapes{%d};',ii));
    eval(sprintf('Phi1=shapes{%d}.Phi;',ii));
    eval(sprintf('A1=shapes{%d}.A;',ii));
    
    for jj=(ii+1):Ns
        eval(sprintf('shape2=shapes{%d};',jj));
        eval(sprintf('Phi2=shapes{%d}.Phi;',jj));
        eval(sprintf('A2=shapes{%d}.A;',jj));
        %%
        eval(sprintf('L120=L{%d,%d};',ii,jj));
        
        %% calculation
%         eval(sprintf('C=X.X%d * X.X%d'';',ii,jj));
        eval(sprintf('C=CC{%d,%d};',ii,jj));
%         C=C.';
        %% point-wise correspondence
%         [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, Phi1, Phi2, 'debug', 0,'numRefinements', 0);
        [shape1ToShape2, shape2ToShape1] = calcP2PFromC(shape1, shape2, C, A1*Phi1, A2*Phi2, 'debug', 0,'numRefinements', 0);
        L12 = [ (1:N)', shape1ToShape2(:) ];
        [ keval_, x ] = kimeval_final2( shape1, shape2, L120, L12(1:fsampl:end,:) );
        keval=keval+keval_;
    
        cnt=cnt+1;
    end
end
%%
keval=keval/cnt;
area = areanum( x, keval );
end