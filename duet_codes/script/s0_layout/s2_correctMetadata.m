clc
clear
close all


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Reads the electrode definition for the cap CY-261.
elec       = ft_read_sens ( 'CY-261.elc' );

% Transforms the ALS definition into RAS definition.
elec.chanpos = [ -1 +1 +1 ] .* elec.chanpos ( :, [ 2 1 3 ] );
elec.elecpos = [ -1 +1 +1 ] .* elec.elecpos ( :, [ 2 1 3 ] );


% Generates an initial layout.
layout     = ft_prepare_layout ( struct ( 'elec', elec ) );

% Loads the corrected layout.
table      = readtable ( 'ant64dry_layout.xlsx' );

% Centers the layout.
table.X    = table.X - mean ( table.X );
table.Y    = table.Y - mean ( table.Y );

% Makes the layout more or less square.
table.X    = table.X / max ( abs ( table.X ) );
table.Y    = table.Y / max ( abs ( table.Y ) );

% Gets the range for the layout and the outline.
lrange     = max ( sqrt ( table.X .^ 2 + table.Y .^ 2 ) );
orange     = max ( sqrt ( layout.outline {1} ( :, 1 ) .^ 2 + layout.outline {1} ( :, 2 ) .^ 2 ) );

% Scales the layout to fit in the outline.
table.X    = 0.95 * table.X * orange / lrange;
table.Y    = 0.95 * table.Y * orange / lrange;

% Imports the corrected coordinates.
hits       = my_matchstr ( layout.label, table.Label );
layout.pos ( hits, 1 ) = table.X;
layout.pos ( hits, 2 ) = table.Y;

ft_plot_layout ( layout, 'outline', 'yes' )


% Generates the neighbour definition.
cfg.method = 'triangulation';
cfg.layout = layout;

neighbours = ft_prepare_neighbours ( cfg );


% Stores the variables.
sensmeta            = [];
sensmeta.system     = 'ANT Dry 64';
sensmeta.label      = elec.label;
sensmeta.elec       = elec;
sensmeta.layout     = layout;
sensmeta.neighbours = neighbours;

% Saves the new layout.
save ( '-v6', 'ant64dry.mat', '-struct', 'sensmeta' )
