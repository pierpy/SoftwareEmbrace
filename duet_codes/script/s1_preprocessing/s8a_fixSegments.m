clc
clear
close all

% Sets the paths.
config.path.segs = '../../data/segments_split/';
config.path.fixs = '../../data/segments_rec/';
config.path.patt = '*.mat';


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions_lap/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Loads the "ideal" electrode definition.
sensmeta = load ( 'ant64dry_hyper' );


% Creates and output folder, if needed.
if ~exist ( config.path.fixs, 'dir' ), mkdir ( config.path.fixs ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.segs, config.path.patt ) );

% Goes through each file.
for findex = 1: numel ( files )
    
    % Pre-loads the data.
    sensdata           = load ( sprintf ( '%s%s', config.path.segs, files ( findex ).name ), 'subject', 'task', 'channel' );
    
    fprintf ( 1, 'Working with subject %s, channel %s.\n', sensdata.subject, sensdata.channel );
    
    
    fprintf ( 1, '  Loading the data.\n' );
    
    % Loads the data.
    sensdata           = load ( sprintf ( '%s%s', config.path.segs, files ( findex ).name ) );
    
    % Gets the trial data.
    trialdata          = sensdata.trialdata;
    
    
    fprintf ( 1, '  Interpolating the missing channels.\n' );
    
    % Gets only the relevant channels (*_1 or *_2).
    rpattern           = sprintf ( '%s$', sensdata.channel ( end - 1: end ) );
    rindex             = ~cellfun ( @isempty, regexp ( sensmeta.label, rpattern ) );
    
    elec               = my_fixsens ( sensmeta.elec, sensmeta.label ( rindex ) );
    
    % Gets the interpolation matrix (spatial filter for the potential).
    chindex            = ismember ( elec.label, trialdata.label );
    wPot               = my_sphsplint ( elec.chanpos ( chindex, : ), elec.chanpos );
    
    
    % Gets the labels for the data (with interpolated channels).
    xlabel             = setdiff ( trialdata.label, elec.label, 'stable' );
    ilabel             = intersect ( trialdata.label, elec.label, 'stable' );
    ixlabel            = trialdata.label;
    olabel             = elec.label;
    oxlabel            = union ( elec.label, trialdata.label, 'stable' );
    
    % Initializes the mapping matrix.
    mapping            = zeros ( numel ( oxlabel ), numel ( ixlabel ) );
    
    % Fills the mapping matrix with the interpolation matrix.
    iindex             = my_matchstr ( ixlabel, ilabel );
    oindex             = my_matchstr ( oxlabel, olabel );
    mapping ( oindex, iindex ) = wPot;
    
    % Fills the mapping matrix with the extended matrix.
    ixindex            = my_matchstr ( ixlabel, xlabel );
    oxindex            = my_matchstr ( oxlabel, xlabel );
    mapping ( oxindex, ixindex ) = eye ( numel ( xlabel ) );
    
%     % Keeps the good channels.
%     igindex            = my_matchstr ( ixlabel, ilabel );
%     ogindex            = my_matchstr ( oxlabel, ilabel );
%     mapping ( ogindex, igindex ) = eye ( numel ( ilabel ) );
    
    
    % Applies the mapping matrix to each trial.
    for tindex = 1: numel ( trialdata.trial )
        trialdata.trial { tindex } = mapping * trialdata.trial { tindex };
    end
    
    % Stores the new channel list.
    trialdata.label    = oxlabel;
    
    % Stores the interpolated value.
    sensdata.trialdata = trialdata;
    
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s', config.path.fixs, files ( findex ).name ), '-struct', 'sensdata' );
end
