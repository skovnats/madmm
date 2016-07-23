function [r] = approxgeodradius(shape)
if isfield( shape, 'idx' )
    idx=shape.idx;
else 
    idx=fps(1,150-1,shape);
end
%%
mx=0;
for i=1:length(idx)
   d = d_shape(shape, idx(i)); 
   mx=max(mx,max(d));
end
r=mx;
end