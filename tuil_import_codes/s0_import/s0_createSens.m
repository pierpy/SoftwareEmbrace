clc
clear
close all

% Defines the system and output file names.
config.system    = 'ANT 256 dry';
config.rawsens   = 'ANT-256ch-dry.elc';
config.file      = 'ANT-256ch-dry.mat';


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% % Loads the raw electrode definitions.
% eegdef = load ( config.rawsens );
% 
% 
% % Initializes the electrode definition.
% % EEGLAB uses ALS and we use RAS.
% elec          = [];
% elec.label    = strtrim ( { eegdef.chanlocs.labels }' );
% elec.chanpos  = cat ( 2, -cat ( 1, eegdef.chanlocs.Y ), cat ( 1, eegdef.chanlocs.X ), cat ( 1, eegdef.chanlocs.Z ) );
% elec.elecpos  = elec.chanpos;
% elec.chantype = repmat ( { 'EEG' }, numel ( elec.label ), 1 );
% elec.chanunit = repmat ( { 'cm' }, numel ( elec.label ), 1 );

% Loads the raw electrode definitions.
elec          = ft_read_sens ( config.rawsens );

% Converts the electrode definition into SI units (meters).
elec          = ft_convert_units ( elec, 'm' );


% Generates a neighbour and layout structures.
cfg                = [];
cfg.elec           = elec;
cfg.method         = 'triangulation';
cfg.projection     = 'polar';

layout             = ft_prepare_layout ( cfg );
neighbours         = ft_prepare_neighbours ( cfg );

% Fixes the layout.
layout.pos ( :, 2 ) = 1.10 * layout.pos ( :, 2 );
layout.pos ( :, 1 ) = 1.10 * layout.pos ( :, 1 );

% Prepares the output.
eegmeta            = [];
eegmeta.system     = config.system;
eegmeta.label      = elec.label;
eegmeta.elec       = elec;
eegmeta.layout     = layout;
eegmeta.neighbours = neighbours;

save ( '-v6', config.file, '-struct', 'eegmeta' );
