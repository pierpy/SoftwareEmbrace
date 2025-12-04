function template = my_sources2roi ( template, roi )


if nargin > 1
    template = my_atlasclean ( template, roi );
end


% Extracts the grid and the atlas.
grid  = template.grid;
atlas = template.atlas;


% Gets the original definitions.
oarea = grid.area ( grid.area > 0 );
oname = atlas.name;
onick = atlas.nick;

% Redefines the grid labeling.
grid.area ( grid.area > 0 ) = 1: sum ( grid.area > 0 );


% Cleans the atlas.
atlas = rmfield ( atlas, setdiff ( fieldnames ( atlas ), { 'atlas' 'name' 'nick' } ) );

% Redefines the area labels.
atlas.name = cellfun ( @(x) sprintf ( 'area%04i', x ), num2cell ( grid.area ( grid.area > 0 ) ), 'UniformOutput', false );
atlas.nick = atlas.name;

% Adds the original labeling.
atlas.oarea = oarea;
atlas.oname = oname;
atlas.onick = onick;


% Stores the modified grid and atlas definitions.
template.grid  = grid;
template.atlas = atlas;
