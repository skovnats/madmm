function  [f] = deltai(n, i)
% creates a delta function of size n with 1 at i

if length( i ) > 1
    f = sparse(n,length(i));
    for j = 1:length(i)
        ii = i(j);
        if ii
            f(ii,j) = 1;
        end
    end
else
    f =zeros(n,1);
    if i
        f(i)=1;
    end
end

end