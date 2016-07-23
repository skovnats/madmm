function [indname,indnamecell] = vec2str(ind)
% Converts vector to a string 4 disp

indname = [];
for i = 1:length( ind )
    if i == 1
        indname = [ indname, [ num2str(ind(i))] ];
        indnamecell{i}=num2str(ind(i));
    else
        indname = [ indname, [ '-' num2str(ind(i))] ];
        indnamecell{i}=num2str(ind(i));
    end
end