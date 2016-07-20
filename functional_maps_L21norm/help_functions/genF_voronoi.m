function [F] = genF_voronoi(shapes, L120, L210, n, show)
% Genrates n corresponding functions on the shapes based on Voronoi
% tesselation

%
N=length( shapes{1}.X );
M=length( shapes{2}.X );

if ~exist( 'show', 'var' )
    show = 1;
end

%-
shape=shapes{1};


% Voronoi Regions
%-fps
% ind=fps(1,n-1,shape); ind=ind(:);
np = randperm(N);
ind=fps(np(1),n-1,shape); ind=ind(:);
ind2=L120(ind,2);

%-
F{1}= voronoiregions2( shape, ind );
% F{1}= voronoiregions( shape, ind );

%-
%%%%%% !!!!! NOTE !!!! %%%%%%%
% This is true only for F{2}
% F{2} = F{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F{2}=zeros(M,n);
for i = 1:n
    ind=find(F{1}(:,i));
    ind=ind(:);
    %
    ind0=ind;
    %
    ind=L120(ind,2);
    
    %
    F{2}(ind,i) = 1;
    
    %
    ind = ismember( L210(:,2), ind0  );
    ind = find(ind);
    F{2}(ind,i) = 1;
    
    %
%     d=d_shape2(shapes{2},ind);
%     ind12=d<2.5;
%     
%     %
%     ind21=ismember( L210(:,2), ind0 );
%     
%     %
%     ind=logical(ind21+ind12);
%     
%     %
%     F{2}(ind,i) = 1;
end

if show
    dshr(shapes{1}, F{1})
    figure
    dshr(shapes{2}, F{2})
end
end