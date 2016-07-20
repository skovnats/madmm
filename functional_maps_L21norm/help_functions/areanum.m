function [ area ] = areanum( x,f )
%numerically calc area/int under the curve

area=0;
for i=1:(length(x)-1)
    area=area + f(i)*(x(i+1)-x(i));
end

end

