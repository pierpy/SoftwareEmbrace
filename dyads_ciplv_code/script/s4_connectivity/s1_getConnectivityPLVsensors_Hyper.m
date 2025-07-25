clc
clear
close all

% Sets the paths.
config.path.sens  = '../../data/segments/';
config.path.conn  = '../../data/connectivity/plv/';
config.path.patt  = '*.mat';

% Action when the task have already been processed.
config.overwrite  = true;

% Defines the bands.
config.bands (1).name   = 'Theta';
config.bands (1).length = 1.8;
config.bands (1).edges  = [  4.0  8.0 ];

config.bands (2).name   = 'Alpha';
config.bands (2).length = 1.8;
config.bands (2).edges  = [  8.0 12.0 ];

config.bands (3).name   = 'Low-beta';
config.bands (3).length = 1.8;
config.bands (3).edges  = [ 12.0 20.0 ];

config.bands (4).name   = 'High-beta';
config.bands (4).length = 1.8;
config.bands (4).edges  = [ 20.0 30.0 ];

config.bands (5).name   = 'Beta';
config.bands (5).length = 1.8;
config.bands (5).edges  = [ 12.0 30.0 ];

config.bands (6).name   = 'Gamma';
config.bands (6).length = 1.8;
config.bands (6).edges  = [ 30.0 45.0 ];

% Defines the parameters.
config.downsample = false;
config.keeptrials = true;


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Creates and output folder, if needed.
if ~exist ( config.path.conn, 'dir' ), mkdir ( config.path.conn ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.sens, config.path.patt ) );

% Goes through each subject.
for findex = 1: numel ( files )
    
    % Pre-loads the data.
    sensdata   = load ( sprintf ( '%s%s', config.path.sens,  files ( findex ).name ), 'subject', 'task', 'stage', 'channel' );
    
    if exist ( sprintf ( '%s%s_%s%s_%s.mat', config.path.conn, sensdata.subject, sensdata.task, sensdata.stage, sensdata.channel ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (already calculated).\n', sensdata.subject, sensdata.task, sensdata.channel );
        continue
    end
    
    fprintf ( 1, 'Working on subject ''%s'', task ''%s'', channel group ''%s''.\n', sensdata.subject, sensdata.task, sensdata.channel );
    
    fprintf ( 1, '  Loading data.\n' );
    
    % Loads the data.
    sensdata   = load ( sprintf ( '%s%s', config.path.sens, files ( findex ).name ) );
    
    % Gets the number of padding samples of the data.
    padding    = round ( sensdata.trialinfo.trialpad (1) * sensdata.trialdata.fsample );
    
    
    % Reserves memory for the band data.
    banddatas  = cell ( numel ( config.bands ), 1 );
    
    % Goes through each band.
    for bindex = 1: numel ( config.bands )
        
        % Gets the band information.
        banddata     = config.bands ( bindex );
        
        fprintf ( 1, '  Working in ''%s'' band (%0.1f to %0.1f Hz).\n', banddata.name, banddata.edges );
        
        
        fprintf ( 1, '    Filtering the data.\n' );
        
        % Gets the trial data.
        trialdata    = sensdata.trialdata;
        
        % If the lower edge is above 1 Hz uses a band-pass filter.
        if banddata.edges (1) >= 1
            
            % Desings a band-pass FIR filter.
            filtlen      = round ( banddata.length * trialdata.fsample );
            fir          = fir1 ( filtlen, banddata.edges / ( trialdata.fsample / 2 ) );
            
            % Apllies the filter.
            trialdata    = myft_filtfilt ( fir, 1, trialdata, true );
            
        % Otherwise combines a high-pass IIR and a low-pass FIR filters.
        else
            
            % Designs a low-pass FIR filter.
            filtlen      = round ( banddata.length * trialdata.fsample ) - 2;
            fir          = fir1 ( filtlen, banddata.edges (2) / ( trialdata.fsample / 2 ) );
            
            % Designs a high-pass IIR filter.
            [ iirb, a ]  = butter ( 2, banddata.edges (1) / ( trialdata.fsample / 2 ), 'high' );
            
            % Combines both filters.
            b            = conv ( iirb, fir );
            
            % Applies the combined filter.
            trialdata    = myft_filtfilt ( b, a, trialdata, true );
        end
        
        % Removes the padding.
        trialdata    = myft_rmpad ( trialdata, sensdata.trialinfo.trialpad );
        
        % Extracts the channel data.
        rawdata      = cat ( 3, trialdata.trial {:} );
        
        
        % For hyperscanning:
        % Removes the mean by subject.
        chindex1     = strncmp ( trialdata.label, 'A_', 2 );
        chindex2     = strncmp ( trialdata.label, 'B_', 2 );
        rawdata ( chindex1, : ) = rawdata ( chindex1, : ) - mean ( rawdata ( chindex1, : ), 1 );
        rawdata ( chindex2, : ) = rawdata ( chindex2, : ) - mean ( rawdata ( chindex2, : ), 1 );
        
        
        % Downsamples the data, if requested.
        if config.downsample
            rawdata      = rawdata ( :, floor ( config.downsample / 2 ): config.downsample: end, : );
        end
        
        % Gets the data dimensions.
        nchans       = size ( rawdata, 1 );
        nsamples     = size ( rawdata, 2 );
        ntrials      = size ( rawdata, 3 );
        
        
        fprintf ( 1, '    Calculating sensor-space connectivity.\n' );
        
        % Underlying logic:
        % The more time-consuming calculation is the complex exponential.
        % This algorithm works in the exponential space the whole time.
        % Uses the logarithmic operations for the phase difference:
        % exp ( 1i * ( a - b ) ) = exp ( 1i * a ) / exp ( 1i * b )
        % Performs all the divisions at once with matrix multiplication.
        
        % Reserves memory for the complex matrix.
        cPLV_all     = complex ( zeros ( nchans, nchans, ntrials, 'single' ) );
        
        % Operates for each trial.
        for tindex = 1: ntrials
            
            % Gets the sources time-series for the current trial.
            sourcedata   = rawdata ( :, :, tindex );
            
            % Gets the vector of phases for each time series.
            sourcevector = sourcedata ./ abs ( sourcedata );
            
            % Calculates the PLV as a matrix multiplication.
            cPLV_all ( :, :, tindex )   = sourcevector * sourcevector' / nsamples;
        end
        
        % Looks for the diagonal of the matrix.
        diags        = diag ( nan ( nchans, 1 ) );
        
        % Sets the diagonal to NaN.
        cPLV_all     = cPLV_all + diags;
        
        
        % Calculates the complex and imaginary parts of PLV.
        plv_sens     = abs ( cPLV_all );
        iplv_sens    = abs ( imag ( cPLV_all ) );
        rplv_sens    = abs ( real ( cPLV_all ) );
        ciplv_sens   = iplv_sens ./ sqrt ( 1 - rplv_sens .^ 2 );
        
        
        % Averages along trials, if required.
        if ~config.keeptrials
            plv_sens     = mean ( plv_sens, 3, 'omitnan' );
            ciplv_sens   = mean ( ciplv_sens, 3, 'omitnan' );
        end
        
        
        % Modifies the band structure.
        banddata.label   = trialdata.label;
        banddata.plv     = plv_sens;
        banddata.ciplv   = ciplv_sens;
        
        % Stores the band structure.
        banddatas { bindex } = banddata;
    end
    
    % Concatenates all the bands.
    banddatas        = cat ( 1, banddatas {:} );
    
    
    fprintf ( 1, '  Saving the connectivity data.\n' );
    
    % Prepares the output.
    conndata           = [];
    conndata.subject   = sensdata.subject;
    conndata.task      = sensdata.task;
    conndata.stage     = sensdata.stage;
    conndata.channel   = sensdata.channel;
    conndata.band      = banddatas;
    conndata.trialinfo = sensdata.trialinfo;
    
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s.mat', config.path.conn, conndata.subject, conndata.task, conndata.stage, conndata.channel ), '-struct', 'conndata' );
end
