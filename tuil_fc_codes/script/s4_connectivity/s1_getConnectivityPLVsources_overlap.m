clc
clear
close all

% Sets the paths.
config.path.meg   = '../../data/segments/';
config.path.beam  = '../../data/sources_overlap/beamformers/';
config.path.conn  = '../../data/connectivity_overlap/plv_sources/';
config.path.patt  = '*.mat';

% Action when the task have already been processed.
config.overwrite  = false;

% Defines the template to use.
config.template   = 'template_AAL';

% Defines the bands to use.
config.bands      = { 'Theta' 'Alpha' 'Low-beta' 'High-beta' 'Beta' 'Gamma' };

% Defines the filter to use (or false to use the band-specific filter).
config.filter     = 'Broadband';

% Defines the parameters.
config.downsample = false;
config.keeptrials = false;
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

% Converts each cortical source in its own area.
template = my_sources2roi ( template );

% Gets the number of sources and areas.
nareas   = numel ( template.atlas.name );
nsources = numel ( template.grid.area ( template.grid.area ~= 0 ) );


% Creates and output folder, if needed.
if ~exist ( config.path.conn, 'dir' ), mkdir ( config.path.conn ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.beam, config.path.patt ) );

% Goes through each subject.
for findex = 1: numel ( files )
    
    % Loads the MEG data and the beamformer information.
    filename   = files ( findex ).name;
    [ ~, basename ] = fileparts ( filename );
    
    beamdata   = load ( sprintf ( '%s%s', config.path.beam,  basename ), 'subject', 'task', 'stage', 'channel', 'whitener', 'lambda' );
    
    if exist ( sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s.mat', config.path.conn, beamdata.subject, beamdata.task, beamdata.stage, beamdata.channel, beamdata.whitener, beamdata.lambda ( 1: end -1 ), template.atlas.atlas ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', stage ''%s'', channel group ''%s'', whitening ''%s'', regularization %s (already calculated).\n', beamdata.subject, beamdata.task, beamdata.stage, beamdata.channel, beamdata.whitener, beamdata.lambda );
        continue
    end
    
    fprintf ( 1, 'Working on subject ''%s'', task ''%s'', stage ''%s'', channel group ''%s'', whitening ''%s'', regularization %s.\n', beamdata.subject, beamdata.task, beamdata.stage, beamdata.channel, beamdata.whitener, beamdata.lambda );
    
    fprintf ( 1, '  Loading data.\n' );
    
    megfile    = dir ( sprintf ( '%s%s_%s%s*.mat', config.path.meg, beamdata.subject, beamdata.task, beamdata.stage ) );
    megfile    = megfile (1).name;
    
    sensdata   = myft_load ( sprintf ( '%s%s', config.path.meg,  megfile ) );
    beamdata   = load ( sprintf ( '%s%s', config.path.beam, basename ) );
    
    % Gets the number of padding samples of the data.
    padding    = round ( sensdata.trialinfo.trialpad (1) * sensdata.trialdata.fsample );
    
    
    % Gets the beamformer filter for all the bands, if requested.
    if config.filter
        
        fprintf ( 1, '  Calculating the number of directions for each source.\n' );
        
        % Extracts the selected beamformer.
        bffilter     = strcmpi ( { beamdata.band.name }, config.filter );
        source       = beamdata.band ( bffilter ).sources;
        
        % Gets the requested orientations for each source.
        source       = my_selectOri ( source, config.threshold );
        
        % Keeps only the sources with a valid area identifier.
        areas        = template.grid.inside & template.grid.area > 0;
        filters      = source.filter ( areas );
        sourcearea   = template.grid.area ( areas );
        sourcepos    = template.grid.pos  ( areas, : );
        
        fprintf ( 1, '    Sources with one direction: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 1 ) );
        fprintf ( 1, '    Sources with two directions: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 2 ) );
        fprintf ( 1, '    Sources with three directions: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 3 ) );
        
        
        % Lists the number of orientations for source position.
        nori         = cellfun ( @(x) size ( x, 1 ), filters );
        
        % Sets the source information in cell form.
        nori         = num2cell ( nori );
        sourcearea   = num2cell ( sourcearea );
        sourcepos    = num2cell ( sourcepos, 2 );
        
        % Repeats the area and position for each source.
        sourcearea   = cellfun ( @(x,y) repmat ( x, y, 1 ), sourcearea, nori, 'UniformOutput', false );
        sourcepos    = cellfun ( @(x,y) repmat ( x, y, 1 ), sourcepos,  nori, 'UniformOutput', false );
        
        % Goes back to vector form.
        sourcearea   = cat ( 1, sourcearea {:} );
        sourcepos    = cat ( 1, sourcepos  {:} );
        
        % Writes the beamformer filter in matrix form.
        filter       = cat ( 1, filters {:} );
    end
    
    
    % Keeps only the requested bands.
    bandinfos  = beamdata.band;
    bandinfos  = bandinfos ( ismember ( { bandinfos.name }, config.bands ) );
    
    % Reserves memory for the band data.
    banddatas  = cell ( numel ( bandinfos ), 1 );
    
    % Goes through each band.
    for bindex = 1: numel ( bandinfos )
        
        % Gets the band information.
        banddata     = bandinfos ( bindex );

        fprintf ( 1, '  Working in ''%s'' band (%0.1f to %0.1f Hz).\n', banddata.name, banddata.edges );
        
        
        % Gets the beamformer filter for the band, if requested.
        if ~config.filter
            
            fprintf ( 1, '    Calculating the number of directions for each source.\n' );
            
            % Gets the beamformer filter.
            source       = bandinfos ( bindex ).sources;
            
            % Gets the requested orientations for each source.
            source       = my_selectOri ( source, config.threshold );
            
            % Keeps only the sources with a valid area identifier.
            areas        = template.grid.inside & template.grid.area > 0;
            filters      = source.filter ( areas );
            sourcearea   = template.grid.area ( areas );
            sourcepos    = template.grid.pos  ( areas, : );
            
            fprintf ( 1, '      Sources with one direction: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 1 ) );
            fprintf ( 1, '      Sources with two directions: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 2 ) );
            fprintf ( 1, '      Sources with three directions: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 3 ) );
            
            
            % Lists the number of orientations for source position.
            nori         = cellfun ( @(x) size ( x, 1 ), filters );
            
            % Sets the source information in cell form.
            nori         = num2cell ( nori );
            sourcearea   = num2cell ( sourcearea );
            sourcepos    = num2cell ( sourcepos, 2 );
            
            % Repeats the area and position for each source.
            sourcearea   = cellfun ( @(x,y) repmat ( x, y, 1 ), sourcearea, nori, 'UniformOutput', false );
            sourcepos    = cellfun ( @(x,y) repmat ( x, y, 1 ), sourcepos,  nori, 'UniformOutput', false );
            
            % Goes back to vector form.
            sourcearea   = cat ( 1, sourcearea {:} );
            sourcepos    = cat ( 1, sourcepos  {:} );
            
            % Writes the beamformer filter in matrix form.
            filter       = cat ( 1, filters {:} );
        end
        
        
        fprintf ( 1, '    Filtering the data.\n' );
        
        % Selects only the channels present in the leadfield.
        cfg          = [];
        cfg.channel  = source.label;
        cfg.trials   = 'all';
        
        trialdata    = ft_selectdata ( cfg, sensdata.trialdata );
        
        
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
        
        
        % Gets the order of the channels in the beamformer filter.
        chanorder    = my_matchstr ( trialdata.label, source.label );
        
        % Extracts and sorts the channel data.
        trialdata     = cat ( 3, trialdata.trial {:} );
        trialdata     = trialdata ( chanorder, :, : );
        
        
        % Downsamples the data, if required.
        if config.downsample
            trialdata     = trialdata ( :, floor ( config.downsample / 2 ): config.downsample: end, : );
        end
        
        % Gets the data dimensions.
        nsources     = size ( filter, 1 );
        nchans       = size ( trialdata, 1 );
        nsamples     = size ( trialdata, 2 );
        ntrials      = size ( trialdata, 3 );
        
        
        fprintf ( 1, '    Calculating source-space connectivity.\n' );
        
        % Underlying logic:
        % The more time-consuming calculation is the complex exponential.
        % This algorithm works in the exponential space the whole time.
        % Uses the logarithmic operations for the phase difference:
        % exp ( 1i * ( a - b ) ) = exp ( 1i * a ) / exp ( 1i * b )
        % Performs all the divisions at once with matrix multiplication.
        
        % Memory reservation.
        cPLV_all     = complex ( zeros ( nsources, nsources, ntrials, 'single' ) );
        
        % Operates for each trial.
        for tindex = 1: ntrials
            
            % Gets the sources time-series for the current trial.
            sourcedata   = filter * trialdata ( :, :, tindex );
            
            % Gets the vector of phases for each time series.
            sourcevector = sourcedata ./ abs ( sourcedata );
            
            % Calculates the PLV as a matrix multiplication.
            cPLV_all ( :, :, tindex )   = sourcevector * sourcevector' / nsamples;
        end
        
        % Looks for the (extended) diagonal of the matrix.
        diags        = cellfun ( @nan, nori, 'UniformOutput', false );
        diags        = blkdiag ( diags {:} );
        
        % Sets the all-sources diagonal to NaN.
        cPLV_all     = cPLV_all + diags;
        
        
        fprintf ( '    Calculating the inter-area connectivity.\n' );
        
        % Calculates the mapping matrix from orientations to areas.
        cases         = histcounts ( sourcearea, 'BinMethod', 'integers' );
        weight        = 1 ./ cases ( sourcearea );
        mapping       = sparse ( sourcearea, 1: nsources, weight );
        
        % Calculates the number of auto-orientations for each area.
        autoori       = mapping * isnan ( diags ) * mapping';
        autoori       = diag ( autoori )';
        
        % Calculates the correction factor for the auto-orientations.
        autocorr      = 1 ./ ( 1 - autoori );
        autocorr      = diag ( autocorr - 1 ) + 1;
        
        
        % Reserves memory for the inter-area PLV matrices.
        plv_rms       = zeros ( nareas, nareas, ntrials, 'single' );
        ciplv_rms     = zeros ( nareas, nareas, ntrials, 'single' );
        
        % Goes through each trial.
        for tindex = 1: ntrials
            
            % Gets the current trial.
            cPLV_t      = cPLV_all ( :, :, tindex );
            
            % Sets the diagonal to zeros.
            cPLV_t ( isnan ( diags ) ) = 0;
            
            % Calculates the real and imaginary part of PLV.
            plv_t       = abs ( cPLV_t );
            iplv_t      = abs ( imag ( cPLV_t ) );
            rplv_t      = abs ( real ( cPLV_t ) );
            
            
            % Calculates the square to get the root-mean-square.
            plv_rms_t   = double ( plv_t .* plv_t );
            iplv_rms_t  = double ( iplv_t .* iplv_t );
            rplv_rms_t  = double ( rplv_t .* rplv_t );
            ciplv_rms_t = iplv_rms_t ./ sqrt ( 1 - rplv_rms_t .^ 2 );
            
            % Averages the activity of each area using the map.
            plv_rms_t   = mapping * plv_rms_t * mapping';
            ciplv_rms_t = mapping * ciplv_rms_t * mapping';
%             iplv_rms_t  = mapping * iplv_rms_t * mapping';
%             rplv_rms_t  = mapping * rplv_rms_t * mapping';
            
            % Corrects the auto-orientations.
            plv_rms_t   = plv_rms_t .* autocorr;
            ciplv_rms_t = ciplv_rms_t .* autocorr;
            
            % Takes the square root to get the root mean square.
            plv_rms_t   = sqrt ( plv_rms_t );
            ciplv_rms_t = sqrt ( ciplv_rms_t );
%             iplv_rms_t  = sqrt ( iplv_rms_t );
%             rplv_rms_t  = sqrt ( rplv_rms_t );
%             ciplv_rms_t = iplv_rms_t ./ sqrt ( 1 - rplv_rms_t .^ 2 );
            
            % Stores the inte-area connectivity for the current trial.
            plv_rms   ( :, :, tindex ) = plv_rms_t;
            ciplv_rms ( :, :, tindex ) = ciplv_rms_t;
        end
        
        clear cPLV_all cPLV_ii cPLV_ij
        
        % Averages along trials, if required.
        if ~config.keeptrials
            plv_rms   = mean ( plv_rms, 3, 'omitnan' );
            ciplv_rms = mean ( ciplv_rms, 3, 'omitnan' );
        end
        
        
        % Modifies the band structure.
        banddata         = rmfield ( banddata, 'sources' );
        banddata         = rmfield ( banddata, 'label' );
        banddata.label   = template.atlas.name;
        banddata.nick    = template.atlas.nick;
        banddata.plv_rms = plv_rms;
        banddata.ciplv_rms = ciplv_rms;
        
        % Stores the band structure.
        banddatas { bindex } = banddata;
    end
    
    % Concatenates all the bands.
    banddatas     = cat ( 1, banddatas {:} );
    
    
    fprintf ( 1, '  Saving the connectivity data.\n' );
    
    % Prepares the output.
    conndata          = [];
    conndata.subject  = beamdata.subject;
    conndata.task     = beamdata.task;
    conndata.stage    = beamdata.stage;
    conndata.channel  = beamdata.channel;
    conndata.whitener = beamdata.whitener;
    conndata.lambda   = beamdata.lambda;
    conndata.atlas    = template.atlas.atlas;
    conndata.band     = banddatas;
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s', config.path.conn, conndata.subject, conndata.task, conndata.stage, conndata.channel, conndata.whitener, conndata.lambda ( 1: end -1 ), conndata.atlas ), '-struct', 'conndata' );
end
