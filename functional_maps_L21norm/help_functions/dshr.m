function [] = dshr( shape, U )
% Function get's shape and the binary matrix U (each column represnts the
% region) and plots it

%-
N = size( U, 2 );
n = ceil( sqrt(N) );

%-
for i = 1:N
   subplot( n, n, i );
   dsh( shape, U( :, i ) );
   %-
   title( num2str(i) );
end