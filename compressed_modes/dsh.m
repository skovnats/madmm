% script for dispalying shape
function [] = dsh( varargin )
% input:
%{
{1} - title
{2} - if save
%}
vector = false;
flag = true;
name = [];
switch nargin
    case 1
        shape = varargin{ 1 };
        
    case 2
        shape = varargin{ 1 };
        tname = varargin{ 2 };
        
        if isnumeric(tname)
            vector = true;
        else
            vector = false;
        end
        
    case 3
        shape = varargin{ 1 };
        tname = varargin{ 2 };
        issave = varargin{ 3 };
        
        if isnumeric(tname)
            vector = true;
        else
            vector = false;
        end
        
        if isstr( issave )
            name = issave;
            issave = false;
        end
        
    case 4
        shape = varargin{ 1 };
        tname = varargin{ 2 };
        vector = true;
        name = varargin{ 3 };
        issave = varargin{ 4 };
end

if iscell( shape )
    flag = false;
    for i = 1:length( shape )
       dsh( shape{ i } );
       title( sprintf( 'shape %d/%d', i, length( shape ) ) );
       waitforbuttonpress;
    end
end


if flag
if ~vector
    % displaying
    if ~isfield( shape, 'C' )
        trisurf( shape.TRIV, shape.X, shape.Y, shape.Z, ones(size((shape.X))) ), ...
    end
else
trisurf( shape.TRIV, shape.X, shape.Y, shape.Z, full(tname) ), ...
end

if isfield( shape, 'C' )
    if ~vector
        try
            %%
            % UPD: 15.11.2011
            lab = [shape.L, shape.a, shape.b];
            lab = colorspace( 'lab->rgb', lab );
            shape.C = lab;
            trisurf( shape.TRIV, shape.X, shape.Y, shape.Z, 1:(length(shape.X)) )
        catch
            
        end
        %%
        colormap(shape.C),
        %         colormap(ones( length(shape.X), 3 )),
    end
    axis off, axis image, shading interp;
    %     lighting phong, camlight('headlight'); % was commented
else
    if ~vector
        colormap(ones( length(shape.X), 3 )),
    end
    axis off, axis image, shading interp, lighting phong, camlight('headlight');
end


switch nargin
    case 2
        if ~vector
            title(tname);
        end
    case 3
        if ~vector
            title(tname);
        end
        %
        if  isstr( name )
            title(name);
        end
        
        if issave
            saveas( gcf, [tname '.png'] )
        end
    case 4
        title( name );
        if issave
            saveas( gcf, [name '.png'] )
        end
end
set(gcf,'Color','w');
cameratoolbar;
%%
% try
% caxis( [-max(abs(tname)), max(abs(tname))] );
% colormap(temp(64));
% catch
% end
%%
end