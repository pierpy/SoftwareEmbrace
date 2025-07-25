clc
% clear
close all

% Sets the paths.                       
config.path.meg  = '../../data/segments/';
config.path.lead = '../../template/leadfield/MARS_CTB-10mm_ANT-256ch-dry.mat';
config.path.filt = '../../data/sources_overlap/beamformers/';
config.path.patt = '*.mat';

% Action when the task has already been processed.
config.overwrite = false;

% Defines the channels to use.
config.channel   = 'EEG';

% Defines the regularization parameter.
config.lambda    = '1%';

% Level of whitening: 0 (no whitening), 1 (scaling), 2 (PCA).
config.whitener  = 1;

% Defines the bands.
config.bands (1).name   = 'Theta';
config.bands (1).length = 0.9;
config.bands (1).edges  = [  4.0  8.0 ];

config.bands (2).name   = 'Alpha';
config.bands (2).length = 0.9;
config.bands (2).edges  = [  8.0 12.0 ];

config.bands (3).name   = 'Low-beta';
config.bands (3).length = 0.9;
config.bands (3).edges  = [ 12.0 20.0 ];

config.bands (4).name   = 'High-beta';
config.bands (4).length = 0.9;
config.bands (4).edges  = [ 20.0 30.0 ];

config.bands (5).name   = 'Beta';
config.bands (5).length = 0.9;
config.bands (5).edges  = [ 12.0 30.0 ];

config.bands (6).name   = 'Gamma';
config.bands (6).length = 0.9;
config.bands (6).edges  = [ 30.0 45.0 ];

config.bands (7).name   = 'Broadband';
config.bands (7).length = 0.9;
config.bands (7).edges  = [  2.0 45.0 ];


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path

% Adds the FreeSrufer toolbox to the path.
ft_hastoolbox ( 'spm8', 1, 1 );


% Creates and output folder, if needed.
if ~exist ( config.path.filt, 'dir' ), mkdir ( config.path.filt ); end

% Defines the whitener type.
switch config.whitener
    case 0, config.wname = 'None';
    case 1, config.wname = 'Scaling';
    case 2, config.wname = 'PCA';
    case 3, config.wname = 'MNE';
end


% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.meg, config.path.patt ) );

% Goes through each subject.
for findex = 1: numel ( files )
    
    % Preloads the data.
    megdata             = load ( sprintf ( '%s%s', config.path.meg,  files ( findex ).name ), 'subject', 'task', 'stage', 'channel' );
    
    if ~config.overwrite && exist ( sprintf ( '%s%s_%s%s_%s_w%s_r%s.mat', config.path.filt, megdata.subject, megdata.task, megdata.stage, megdata.channel, config.wname, config.lambda ( 1: end - 1 ) ), 'file' )
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (already calculated).\n', megdata.subject, megdata.task, megdata.channel );
        continue
    end
    if exist ( config.path.lead, 'file' ) ~= 2 && ~exist ( sprintf ( '%s%s.mat', config.path.lead, megdata.subject ), 'file' )
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (no lead field).\n', megdata.subject, megdata.task, megdata.channel );
        continue
    end
    
    fprintf ( 1, 'Calculating beamformers for subject ''%s'', task ''%s'', channel group ''%s''.\n', megdata.subject, megdata.task, megdata.channel );
    
    
    fprintf ( 1, '  Loading the data.\n' );
    
    % Loads the data.
    megdata             = myft_load ( sprintf ( '%s%s', config.path.meg,  files ( findex ).name ) );
    
    % Loads the individual lead field, if provided, or the template one.
    if exist ( config.path.lead, 'dir' )
        leaddata            = load ( sprintf ( '%s%s.mat', config.path.lead, megdata.subject ) );
    else
        leaddata            = load ( config.path.lead );
    end
    
    
    % Loads the two sets of bad channels.
    wetfile             = sprintf ( '%s_%s%s_%s.mat', megdata.subject, strrep ( megdata.task, 'D-', 'W-' ), megdata.stage, megdata.channel );
    wetmeta             = load ( sprintf ( '%s%s', config.path.meg,  wetfile ), 'chaninfo' );
    dryfile             = sprintf ( '%s_%s%s_%s.mat', megdata.subject, strrep ( megdata.task, 'W-', 'D-' ), megdata.stage, megdata.channel );
    drymeta             = load ( sprintf ( '%s%s', config.path.meg,  dryfile ), 'chaninfo' );
    badchan             = union ( wetmeta.chaninfo.bad, drymeta.chaninfo.bad );
    
    % Sets the bad channels for this data.
    megdata.chaninfo.bad = badchan;
    
    
    % Initializes the band information.
    banddatas           = cell ( numel ( config.bands ), 1 );
    
    % Goes through each band.
    for band =  1: numel ( config.bands )
        
        % Gets the current band information.
        banddata            = config.bands ( band );
        
        
        fprintf ( 1, '  Working in ''%s'' band (%0.1f to %0.1f Hz).\n', banddata.name, banddata.edges );
        
        % Gets the data and the lead field definition.
        trialdata           = megdata.trialdata;
        grid                = leaddata.grid;
        
        % Gets only the channels both in the data and the lead field.
        channel             = trialdata.label;
        channel             = setdiff ( channel, megdata.chaninfo.bad );
        channel             = intersect ( channel, grid.label, 'stable' );
        channel             = ft_channelselection ( config.channel, channel );
        
        cfg                 = [];
        cfg.channel         = channel;
        
        trialdata           = ft_selectdata ( cfg, trialdata );
        
        
        fprintf ( 1, '    Filtering the data.\n' );
        
        % If the lower edge is above 1 Hz uses a band-pass filter.
        if banddata.edges (1) >= 1
            
            % Desings a band-pass FIR filter.
            filtlen             = round ( banddata.length * trialdata.fsample );
            fir                 = fir1 ( filtlen, banddata.edges / ( trialdata.fsample / 2 ) );
            
            % Apllies the filter.
            trialdata           = myft_filtfilt ( fir, 1, trialdata );
            
        % Otherwise combines a high-pass IIR and a low-pass FIR filters.
        else
            
            % Designs a low-pass FIR filter.
            filtlen             = round ( banddata.length * trialdata.fsample ) - 2;
            fir                 = fir1 ( filtlen, banddata.edges (2) / ( trialdata.fsample / 2 ) );
            
            % Designs a high-pass IIR filter.
            [ iirb, a ]         = butter ( 2, banddata.edges (1) / ( trialdata.fsample / 2 ), 'high' );
            
            % Combines both filters.
            b                   = conv ( iirb, fir );
            
            % Applies the combined filter.
            trialdata           = myft_filtfilt ( b, a, trialdata );
        end
        
        % Removes the padding.
        trialdata           = myft_rmpad ( trialdata, megdata.trialinfo.trialpad );
        
        
        fprintf ( 1, '    Calculating the evoked response and covariance matrix.\n' );
        
        % Calculates the evoked response.
        cfg                 = [];
        cfg.channel         = 'all';
        cfg.trials          = 'all';
        cfg.keeptrials      = false;
        cfg.latency         = [ -Inf +Inf ];
        cfg.baseline        = [ -Inf +Inf ];
        
        timelock            = myft_timelock ( cfg, trialdata );
        
        
%         % Gets the optimal sampling rate for the evoked response.
%         downsamp            = floor ( ( trialdata.fsample / 2 ) / banddata.edges (2) );
%         
%         % Gets the index to the closest-to-zero sample.
%         [ ~, zerosamp ]     = min ( abs ( timelock.time ) );
%         firstsamp           = rem ( zerosamp - 1, downsamp ) + 1;
%         
%         % Downsamples the evoked response.
%         timelock.avg        = timelock.avg  ( :, firstsamp: downsamp: end );
%         timelock.var        = timelock.var  ( :, firstsamp: downsamp: end );
%         timelock.time       = timelock.time ( :, firstsamp: downsamp: end );
%         timelock.dof        = timelock.dof  ( :, firstsamp: downsamp: end );
        
        
        fprintf ( 1, '    Creating the beamformer spatial filter.\n' );
        
        % Generates the data whitener, if requested.
        whitener            = my_whitener ( timelock, config.whitener, true );
        
        % Calculates the spatial filter.
        cfg                 = [];
        cfg.grid            = grid;
        cfg.keepori         = 'yes';
        cfg.keepleadfield   = 'no';
        cfg.method          = 'lcmv';
        
        cfg.lcmv.lambda     = config.lambda;
        cfg.lcmv.projectmom = 'no';
        cfg.lcmv.keepmom    = 'no';
        cfg.lcmv.keepcov    = 'yes';
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.feedback   = 'no';
        cfg.lcmv.subspace   = whitener;
        
        sources             = my_sourceanalysis ( cfg, timelock );
        
        
        % Sets the output
        banddata.label      = sources.label;
        banddata.sources    = sources;
        
        banddatas { band }  = banddata;
    end
    
    % Rewrites the band informations as an array of structures.
    banddatas           = cat ( 1, banddatas {:} );
    
    
    fprintf ( 1, '  Saving the data.\n' );
    
    % Prepares the output.
    sourcedata          = [];
    sourcedata.subject  = megdata.subject;
    sourcedata.task     = megdata.task;
    sourcedata.stage    = megdata.stage;
    sourcedata.channel  = config.channel;
    sourcedata.whitener = config.wname;
    sourcedata.lambda   = strrep ( config.lambda, '%', 'pc' );
    sourcedata.band     = banddatas;
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s_w%s_r%s', config.path.filt, sourcedata.subject, sourcedata.task, sourcedata.stage, sourcedata.channel, sourcedata.whitener, sourcedata.lambda ), '-struct', 'sourcedata' );
end
