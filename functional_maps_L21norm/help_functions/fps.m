function [sample, d] = fps(fstSampl, Ns, shape,s)
% D - is a distance matrix between sampled points.

if nargin < 4
%     s=distance_calc('init', 'geodesic',shape);
    s=dstcalc('init', 'geodesic',shape);
end
% s=distance_calc('init', 'euc',shape);

Nv = s.Nv;
sample = [fstSampl];
initnumofsamples = length(sample);
% ds = zeros(Nv, 1) + Inf;
% ds = distance_calc('compute', s, sample,shape);
ds = dstcalc('compute', s, sample,shape);
str = '';
for k=1:Ns,
    
    % Compute distance map from current sources
%     d = distance_calc('compute', s, sample(k-1+initnumofsamples),shape);
    d = dstcalc('compute', s, sample(k-1+initnumofsamples),shape);
    ds = min(ds, d);
    str = mprintf(str, 'Sampling %d / %d', k, Ns);
    
%         if nargin > 2,
%             trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, ds); axis image;
%             shading interp;
%             lighting phong;
%             camlight head;
%             axis off;
%             hold on;
%             plot3(shape.X(sample), shape.Y(sample), shape.Z(sample), '.r');
%             hold off;
%             drawnow;
%         end
    
    [r,sample(k+initnumofsamples)] = max(ds);
end
%-
sample = sample(1:(Ns+initnumofsamples));


%{
for k=1:Ns,
    
    % Compute distance map from current sources
    d = distance_calc('compute', s, sample(k),shape);
    D(:,k) = d(sample);
    str = mprintf(str, 'Distance map %d / %d', k, Ns);
    
    if nargin > 1,
    trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, d); axis image;
    shading interp;
    lighting phong;
    camlight head;
    axis off;
    hold on;
    plot3(shape.X(sample(k)), shape.Y(sample(k)), shape.Z(sample(k)), '.r');
    hold off;
    drawnow;
    end
    
end
%}
str = mprintf(str, '');

