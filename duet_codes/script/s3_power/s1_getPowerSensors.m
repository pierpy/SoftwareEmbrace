clc
clear
close all

% Paths IN/OUT.
config.path.meg   = '../../data/segments/';
config.path.pow   = '../../data/spectra/dpss_05_new/';
config.path.patt  = '*.mat';

% Frequency band to calculate.
config.band       = [  2 45 ];

% Taper to use.
config.taper      = 'dpss';
config.tapsmofrq  = 0.5;

% Sets the action when the task has already been processed.
config.overwrite  = false;


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'spm8', 1, 1 );


% Generates the output folder, if needed.
if ~exist ( config.path.pow, 'dir' ), mkdir ( config.path.pow ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.meg, config.path.patt ) );

% Goes through each subject.
for file = 1: numel ( files )
    
    % Pre-loads the data.
    sensdata         = load ( sprintf ( '%s%s', config.path.meg, files ( file ).name ), 'subject', 'task', 'stage', 'channel' );
    
    % Gets the text of the message.
    msgtext          = sprintf ( 'subject ''%s'', task ''%s''', sensdata.subject, sensdata.task );
    if ~isempty ( sensdata.stage )
        msgtext        = sprintf ( '%s, stage ''%s''', msgtext, sensdata.stage );
    end
    msgtext          = sprintf ( '%s, channel group ''%s''', msgtext, sensdata.channel );
    
    if exist ( sprintf ( '%s%s_%s%s_%s.mat', config.path.pow, sensdata.subject, sensdata.task, sensdata.stage, sensdata.channel ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring %s (Already calculated).\n', msgtext );
        continue
    end
    
    fprintf ( 1, 'Working with %s.\n', msgtext );
    
    
    fprintf ( 1, '  Loading the data.\n' );
    
    % Loads the data.
    sensdata         = load ( sprintf ( '%s%s', config.path.meg, files ( file ).name ) );
    trialdata        = sensdata.trialdata;
    
    % Removes the padding.
    trialdata        = myft_rmpad ( trialdata, sensdata.trialinfo.trialpad );
    
    
    fprintf ( 1, '  Calculating the power spectra.\n' );
    
    % Gets the spectra.
    cfg              = [];
    cfg.foilim       = config.band;
    cfg.method       = 'mtmfft';
    cfg.taper        = config.taper;
    cfg.tapsmofrq    = config.tapsmofrq;
    cfg.output       = 'power';
    
    freqdata         = myft_freqanalysis ( cfg, trialdata );
    
    
%     % Gets the spectra.
%     cfg              = [];
%     cfg.channel      = 'all';
%     cfg.trials       = 'all';
%     cfg.foilim       = config.band;
%     cfg.method       = 'mtmfft';
%     cfg.output       = 'power';
%     cfg.taper        = 'dpss';
%     cfg.tapsmofrq    = 2;
%     cfg.keeptrials   = 'no';
%     cfg.feedback     = 'none';
%     
%     freqdata         = ft_freqanalysis ( cfg, trialdata );
    
    
    fprintf ( 1, '  Saving the spectrum.\n' );
    
    % Saves the data.
    powdata          = [];
    powdata.subject  = sensdata.subject;
    powdata.task     = sensdata.task;
    powdata.stage    = sensdata.stage;
    powdata.channel  = sensdata.channel;
    powdata.fileinfo = sensdata.fileinfo;
    powdata.freqdata = freqdata;
    
    save ( '-v6', sprintf ( '%s%s_%s%s_%s', config.path.pow, powdata.subject, powdata.task, powdata.stage, powdata.channel ), '-struct', 'powdata' );
end
