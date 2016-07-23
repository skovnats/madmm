function [ind_rsamples] = rsamples( shape, ind_src, r )
% Function gets shape, index of vertecis, and returns a new random samples withing
% radius r

%- tolerance of the circle
tol = 0.05; % 5%

ind_rsamples = [];
for i = 1:length(ind_src)
    %- current point
    ii = ind_src( i );
    %- dist
    [d] = d_shape(shape, ii);
    %- take
    try
        ind = find( (d <= (r + r * tol)).*(d >= (r - r * tol)) );
        %-
        ind = ind(randperm( length(ind) ));
        %-
        cnt = 1;
        while ~isempty( intersect(ind_rsamples, ind(cnt)) )
            cnt = cnt + 1;
        end
        %-
        ind_rsamples = [ ind_rsamples; ind(cnt) ];
    catch
        ind = find( (d <= (r + r * 2*tol)).*(d >= (r - r * 2*tol)) );
        %-
        ind = ind(randperm( length(ind) ));
        %-
        cnt = 1;
        while ~isempty( intersect(ind_rsamples, ind(cnt)) )
            cnt = cnt + 1;
        end
        %-
        ind_rsamples = [ ind_rsamples; ind(cnt) ];
    end
end