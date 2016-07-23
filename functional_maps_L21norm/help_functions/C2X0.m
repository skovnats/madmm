function [X0]=C2X0(CC,type)

if ~exist('type','var')
    type=1;
end

switch type
    case 1
        C=CC{1,2};
        [u1,d,u2]=svd(C);
        X0.X1=u1;
        X0.X2=u2;
        
        % N=numel(fields(X));
        N=size(CC,1);
        
        for i = 3:N
            [u,~,v] = svd(CC{1,i});
            eval(sprintf('X0.X%d = v * u'' * u1;', i));
        end
    case 2
        N=size(CC,1);
        
        for i = 1:N
            eval(sprintf('X0.X%d = CC{%d,1};', i,i));
        end
end

end
