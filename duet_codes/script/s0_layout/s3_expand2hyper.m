clc
clear
close all


% Loads the original Duke layout.
sensmeta = load ( 'ant64dry.mat' );


% Extracts the electrode definition.
elec       = sensmeta.elec;

% Duplicates the electrode definitions.
elec.chanpos  = cat ( 1, ...
    elec.chanpos, ...
    elec.chanpos + 10 * [ 1 1 1 ] );
elec.chantype = cat ( 1, ...
    elec.chantype, ...
    elec.chantype );
elec.chanunit = cat ( 1, ...
    elec.chanunit, ...
    elec.chanunit );
elec.elecpos  = cat ( 1, ...
    elec.elecpos, ...
    elec.elecpos );

% Duplicates the labels.
elec.label    = cat ( 1, ...
    strcat ( elec.label, '_1' ), ...
    strcat ( elec.label, '_2' ) );


% Extracts the layout.
layout     = sensmeta.layout;

% Duplicates the electrodes in the layout.
layout.pos    = cat ( 1, ...
    layout.pos ( 1: end - 2, : ), ...
    layout.pos ( 1: end - 2, : ), ...
    layout.pos ( end - 1: end, : ) );
layout.width  = cat ( 1, ...
    layout.width ( 1: end - 2 ), ...
    layout.width ( 1: end - 2 ), ...
    layout.width ( end - 1: end ) );
layout.height = cat ( 1, ...
    layout.height ( 1: end - 2 ), ...
    layout.height ( 1: end - 2 ), ...
    layout.height ( end - 1: end ) );

% Duplicates the labels.
layout.label  = cat ( 1, ...
    strcat ( layout.label ( 1: end - 2 ), '_1' ), ...
    strcat ( layout.label ( 1: end - 2 ), '_2' ), ...
    layout.label ( end - 1: end ) );


% Extracts the neighbour definition.
neighbours = sensmeta.neighbours;

% Duplicates the neighbours.
neigh1     = neighbours (:);
neigh2     = neighbours (:);

% Updates the labels of the first set.
for nindex = 1: numel ( neigh1 )
    neigh1 ( nindex ).label       = strcat ( neigh1 ( nindex ).label, '_1' );
    neigh1 ( nindex ).neighblabel = strcat ( neigh1 ( nindex ).neighblabel, '_1' );
end

% Updates the labels of the second set.
for nindex = 1: numel ( neigh2 )
    neigh2 ( nindex ).label       = strcat ( neigh2 ( nindex ).label, '_2' );
    neigh2 ( nindex ).neighblabel = strcat ( neigh2 ( nindex ).neighblabel, '_2' );
end

% Concatenates both sets.
neighbours = cat ( 1, neigh1, neigh2 );


% Stores the variables.
sensmeta            = [];
sensmeta.system     = 'ANT Dry 64 with hyperscanning';
sensmeta.label      = elec.label;
sensmeta.elec       = elec;
sensmeta.layout     = layout;
sensmeta.neighbours = neighbours;

% Saves the new layout.
save ( '-v6', 'ant64dry_hyper.mat', '-struct', 'sensmeta' )
