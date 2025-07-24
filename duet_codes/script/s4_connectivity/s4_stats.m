clc
clear
close all

% Sets the paths.
config.path.conn  = '../../data/connectivity/plv/';
config.path.stats = '../../stats/connectivity/plv_sessions/';
config.path.patt  = '*.mat';

% Defines the metric.
config.metric     = 'ciplv';

% Defines the tasks to compare.
config.comp (1).label = 'playing-(1)';
config.comp (1).task  = 'playing1';
config.comp (2).label = 'resting-(1)';
config.comp (2).task  = 'resting-EO1';
config.comp (3).label = 'playing-(2)';
config.comp (3).task  = 'playing2';
config.comp (4).label = 'resting-(2)';
config.comp (4).task  = 'resting-EO2';


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Creates and output folder, if required.
if ~exist ( config.path.stats, 'dir' ), mkdir ( config.path.stats ); end


% Goes through each comparison.
for cindex = 1: numel ( config.comp )
    
    % Gets the metadata for the current comparison.
    compinfo  = config.comp ( cindex );
    
    fprintf ( 1, 'Performing comparison %s (%s).\n', compinfo.label, compinfo.task )
    
    
    fprintf ( 1, '  Loading the connectivity values.\n' )
    
    % Loads the connectivity values.
    connfile  = dir ( sprintf ( '%s*_%s_*', config.path.conn, compinfo.task ) );
    conndata  = load ( sprintf ( '%s%s', config.path.conn, connfile (1).name ) );
    
    % Goes throug each band.
    for bindex = 1: numel ( conndata.band )
        
        % Gets the metadata for the current band.
        bandinfo  = conndata.band ( bindex );
        
        fprintf ( 1, '    Working on band %s (%.1f-%.1f Hz).\n', bandinfo.name, bandinfo.edges );
        
        
        % Gets the list of trials per each file.
        trials1   = conndata.trialinfo.trialfile == 1;
        trials2   = conndata.trialinfo.trialfile == 2;
        
        % Gets the band data for each task.
        banddata1 = conndata.band ( bindex ).( config.metric ) ( :, :, trials1 );
        banddata2 = conndata.band ( bindex ).( config.metric ) ( :, :, trials2 );
        
        % Calculates the average connectivity for each task.
        connmean1 = nanmean ( banddata1, 3 );
        connmean2 = nanmean ( banddata2, 3 );
        
        % Permutes the data into a shape acceptable for t-test.
        banddata1 = permute ( banddata1, [ 3 1 2 ] );
        banddata2 = permute ( banddata2, [ 3 1 2 ] );
        
        % Performs the comparison.
        [ ~, p, ~, s ] = ttest2 ( banddata1, banddata2 );
        pvalue    = permute ( p, [ 2 3 1 ] );
        tstat     = permute ( s.tstat, [ 2 3 1 ] );
        
        
        % Prepares the output.
        banddata  = [];
        banddata.subject    = conndata.subject;
        banddata.comparison = compinfo.label;
        banddata.task       = compinfo.task;
        banddata.bandname   = bandinfo.name;
        banddata.bandedges  = bandinfo.edges;
        banddata.metric     = config.metric;
        banddata.label      = bandinfo.label;
        banddata.connmean1  = connmean1;
        banddata.connmean2  = connmean2;
        banddata.conndiff   = connmean1 - connmean2;
        banddata.tstat      = tstat;
        banddata.pvalue     = pvalue;
        
        % Saves the data.
        save ( '-v6', sprintf ( '%s%s_%s_%s_%s.mat', config.path.stats, banddata.subject, banddata.comparison, banddata.metric, banddata.bandname ), '-struct', 'banddata' )
    end
end
