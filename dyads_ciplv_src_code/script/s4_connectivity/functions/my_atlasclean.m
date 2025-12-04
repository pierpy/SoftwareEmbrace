function template = my_atlasclean ( template, roi )


if nargin < 2
    roi = unique ( template.grid.area );
    roi = setdiff ( roi, 0 );
end


% Extracts the grid and the atlas.
grid  = template.grid;
atlas = template.atlas;


% Cleans the grid.
grid.area ( ~grid.inside ) = 0;
grid.inside = grid.area > 0;

% Lists the areas both selected and present in the grid.
area  = intersect ( roi, grid.area );

% Cleans the atlas.
atlas = rmfield ( atlas, setdiff ( fieldnames ( atlas ), { 'atlas' 'name' 'nick' 'order' 'pos' 'unit' } ) );
atlas.name  = atlas.name  ( area );
atlas.nick  = atlas.nick  ( area );
atlas.pos   = atlas.pos   ( area, : );

% Re-orders the areas in the atlas.
if isfield ( atlas, 'order' )
    
    % Cleans the atlas.
    atlas.order = atlas.order ( area, : );
    
    % Re-orders the areas.
    [ ~, idx ] = sort ( atlas.order );
    order  = zeros ( size ( atlas.order ) );
    order ( idx ) = 1: numel ( atlas.order );
    
    % Stores the modified order.
    atlas.order = order;
end


% Gets the original grid labeling.
oarea = grid.area;

% Clears the grid labeling.
grid.area (:) = 0;

% Relabels the grid, if required.
for aindex = 1: numel ( area )
    
    % Lists the sources in the current area.
    asource = oarea == area ( aindex );
    
    % Relabels the sources.
    grid.area ( asource ) = aindex;
end


% Stores the modified grid and atlas definitions.
template.grid  = grid;
template.atlas = atlas;
