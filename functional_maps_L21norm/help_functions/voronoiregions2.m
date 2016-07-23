function [ vrf, f ] = voronoiregions2( shape, ind )
%- 
nv = length( shape.X );
vr = zeros( nv, length(ind) );
f = zeros( nv, 1 );

%-
for i = 1:length(ind)
    %
    d = d_shape2( shape, ind( i ) );
    
    %
    vr( :, i ) = d;
end

%
[~,ind_] = min(vr');
vrf = vr;

%
for i = 1:length( ind )
    vrf(:,i) = double( ind_(:) == i ); 
end
end

