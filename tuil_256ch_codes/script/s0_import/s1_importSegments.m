clc
clear
close all

% Defines the paths.
config.path.raw  = '../../data/eeg/';
config.path.seg  = '../../data/segments/';
config.path.patt = '*.cnt';

% Sets the segmentation options.
config.trialfun  = 'restingSegmentation';
config.segment   = 2.0;
config.padding   = 1.0;
config.overlap   = 1.0;


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions_eep', pwd ) );

% Adds FieldTrip to the path, if required.
myft_path


% Creates the output folder, if requried.
if ~exist ( config.path.seg, 'dir' ), mkdir ( config.path.seg ), end


% Lists the files in the raw data folder.
files = dir ( sprintf ( '%s%s', config.path.raw, config.path.patt ) );

% Goes through each file.
for findex = 1: numel ( files )
    
    % Gets the subject and the task.
    dummy       = regexp ( files ( findex ).name, 'EEG([\d])+_(D|W)_(alpha|resting)', 'tokens' );
    subject     = sprintf ( 'EEG%02s', dummy {1} {1} );
    task        = sprintf ( '%s-%s', dummy {1} { 2: 3 } );
    
    fprintf ( 1, 'Working with subject %s, task %s.\n', subject, task );
    
    
    
    % Gets the full file name.
    dataset     = sprintf ( '%s%s', config.path.raw, files ( findex ).name );
    
    
    % Reads the header.
    header              = my_read_header ( dataset );
    
    % Sets the file information.
    fileinfo            = [];
    fileinfo.dataset    = dataset;
    fileinfo.subject    = subject;
    fileinfo.task       = task;
    fileinfo.stage      = '';
    fileinfo.index      = 1;
    fileinfo.begtime    = NaN;
    fileinfo.endtime    = NaN;
    fileinfo.header     = header;
    fileinfo.headshape  = [];
    fileinfo.event      = [];
    
    
    fprintf ( 1, '  Reading the data from disk.\n' );
    
    % Reads the data.
    cfg                 = [];
    cfg.dataset         = dataset;
    cfg.header          = header;
    
    wholedata           = my_read_data ( cfg );
    
    
    fprintf ( 1, '  Segmenting the data into %.2f s segments with %.2f s overlap.\n', config.segment, config.overlap );
    
    % Extracts the trials.
    trialfun              = str2func ( config.trialfun );
    
    fileconfig          = config;
    fileconfig.dataset  = fileinfo.dataset;
    fileconfig.header   = fileinfo.header;
    fileconfig.event    = fileinfo.event;
    fileconfig.begtime  = fileinfo.begtime;
    fileconfig.endtime  = fileinfo.endtime;
    fileconfig.addpadd  = true;
    fileconfig.feedback = 'no';
    
    trialdef            = trialfun ( fileconfig );
    
    % Segments the data.
    fileconfig.feedback = 'no';
    fileconfig.trl      = trialdef;
    
    trialdata           = ft_redefinetrial ( fileconfig, wholedata );
    clear wholedata
    
    
    % Sets the trial information.
    trialinfo.trialdef  = trialdef;
    trialinfo.trialfile = ones ( size ( trialdef, 1 ), 1 );
    trialinfo.trialpad  = config.padding * ones ( size ( trialdef, 1 ), 1 );
    
    
    % Reads the bad channel information.
    metadata    = load ( regexprep ( dataset, '.cnt$', '.mat' ) );
    badchan     = setdiff ( trialdata.label, trialdata.label ( metadata.ch_subset_sorted ) );
    
    % Sets the bad channel information.
    chaninfo            = [];
    chaninfo.bad        = badchan;
    
    
    fprintf ( 1, '  Saving the data.\n' );
    
    % Prepares the epoch data.
    epochdata           = [];
    epochdata.subject   = subject;
    epochdata.task      = task;
    epochdata.stage     = '';
    epochdata.channel   = 'EEG';
    epochdata.fileinfo  = fileinfo;
    epochdata.chaninfo  = chaninfo;
    epochdata.trialinfo = trialinfo;
    epochdata.trialdata = trialdata;
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s.mat', config.path.seg, epochdata.subject, epochdata.task, epochdata.stage, epochdata.channel ), '-struct', 'epochdata' )
end
