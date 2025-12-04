clc
clear
close all

% Sets the paths.
config.path.meg   = '../../data/segments_window/';
config.path.inv   = '../../data/sources/lcmv/';
config.path.conn  = '../../data/connectivity/plv_areas_window/';
config.path.patt  = '*.mat';

% Action when the task have already been processed.
config.overwrite  = false;

% Defines the template to use.
config.template   = 'template_AAL_hyper';
config.atlas      = 'AAL_hyper';

% Defines the bands to use.
config.bands      = { 'Theta' 'Alpha' 'Low-beta' 'High-beta' 'Beta' 'Gamma' };

% Defines the filter to use (or false to use the band-specific filter).
config.filter     = 'Broadband';

% Defines the parameters.
config.downsample = false;
config.keeptrials = true;
config.threshold  = .99;


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'spm8', 1, 1 );


% Loads the template for the label information.
template = load ( config.template, 'grid', 'atlas' );

% Gets the number of sources and areas.
nareas   = numel ( template.atlas.name );
nsources = numel ( template.grid.area ( template.grid.area ~= 0 ) );


% Creates and output folder, if needed.
if ~exist ( config.path.conn, 'dir' ), mkdir ( config.path.conn ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.inv, config.path.patt ) );

% Goes through each subject.
for findex = 1: numel ( files )
    
    % Loads the MEG data and the beamformer information.
    filename   = files ( findex ).name;
    [ ~, basename ] = fileparts ( filename );
    
    invdata    = load ( sprintf ( '%s%s', config.path.inv,  basename ), 'subject', 'task', 'stage', 'channel', 'whitener', 'lambda' );
    
    % if exist ( sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s.mat', config.path.conn, invdata.subject, invdata.task, invdata.stage, invdata.channel, invdata.whitener, invdata.lambda ( 1: end -1 ), template.atlas.atlas ), 'file' ) && ~config.overwrite
    %     fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', stage ''%s'', channel group ''%s'', whitening ''%s'', regularization %s (already calculated).\n', invdata.subject, invdata.task, invdata.stage, invdata.channel, invdata.whitener, invdata.lambda );
    %     continue
    % end
    
    fprintf ( 1, 'Working on subject ''%s'', task ''%s'', stage ''%s'', channel group ''%s'', whitening ''%s'', regularization %s.\n', invdata.subject, invdata.task, invdata.stage, invdata.channel, invdata.whitener, invdata.lambda );
    
    fprintf ( 1, '  Loading data.\n' );
    
    sensfile   = dir ( sprintf ( '%s%s_%s%s*.mat', config.path.meg, invdata.subject, invdata.task, invdata.stage ) );
    sensfile   = sensfile (1).name;
    
    sensdata   = load ( sprintf ( '%s%s', config.path.meg,  sensfile ) );
    invdata    = load ( sprintf ( '%s%s', config.path.inv, basename ) );
    
    conndata   = load ( sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s', config.path.conn, invdata.subject, invdata.task, invdata.stage, invdata.channel, invdata.whitener, invdata.lambda, config.atlas ) );


    fprintf ( 1, '  Adding the trial information.\n' );

    % Adds the trial data to the connectivity data.
    conndata.trialinfo = sensdata.trialinfo;


    fprintf ( 1, '  Saving the data.\n' )
    
    % Saves the data.
    save ( '-v7.3', sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s', config.path.conn, conndata.subject, conndata.task, conndata.stage, conndata.channel, conndata.whitener, conndata.lambda, conndata.atlas ), '-struct', 'conndata' );
end
