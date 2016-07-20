function [CC]=C2C(CC0)
% CC=cell;
for i=1:size(CC0,1)
    for j=1:size(CC0,1)
        CC{i,j} = CC0{1,j} * CC0{i,1};
    end
end
end