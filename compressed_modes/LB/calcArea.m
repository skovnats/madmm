function [area_vert, area_tri] = calcArea(shape)
% Calculate area of each triangle
% area_tri = cross(B - A, C - A, 2);
% area_tri = 1/2*sqrt(sum(area_tri.^2, 2));

Vertices = [shape.X, shape.Y, shape.Z];
A = Vertices(shape.TRIV(:,1), :);
B = Vertices(shape.TRIV(:,2), :);
C = Vertices(shape.TRIV(:,3), :);
area_tri = 1/2*sqrt(sum((B - A).^2, 2).*sum((C - A).^2, 2) - dot(B - A, C - A, 2).^2);

if (nargout == 2)
    N = numel(shape.X);
    area_vert = zeros(N, 1);
    for k = 1:N
        [rows, cols] = find(shape.TRIV == k); %#ok<*NASGU>
        area_vert(k) = sum(area_tri(rows));
    end
else
    area_vert = area_tri;
end

end