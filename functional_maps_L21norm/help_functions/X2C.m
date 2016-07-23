function CC = X2C(X)

N=numel( fieldnames(X) );

for i=1:N
    for j=1:N
        
%         % original
%         eval(sprintf('CC{%d,%d} = X.X%d * X.X%d.'';', i, j, i, j ));
        
        % fixed by E.Rodola, 26.06.2016
        eval(sprintf('CC{%d,%d} = (X.X1*X.X%d'')*(X.X1*X.X%d'')'';', i, j, j, i ));
        
%         figure
%         subplot(121),imagesc(CC{3,5})
%         subplot(122),imagesc((X0.X1*X0.X5')*(X0.X1*X0.X3')')
    end
end

end
