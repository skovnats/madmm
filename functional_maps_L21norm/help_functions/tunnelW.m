function [W] = tunnelW(k)

params.k=k;
W = zeros(params.k);
    for i=1:params.k
        for j=1:params.k
            slope = [1 1]./sqrt(2);
            W(i,j) = exp(-0.03*sqrt(i.^2 + j.^2))*norm(cross([slope 0], [i,j, 0]-[1 1 0]));        
        end
    end
end