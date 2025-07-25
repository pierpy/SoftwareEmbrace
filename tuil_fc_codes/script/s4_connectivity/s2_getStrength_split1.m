clc
clear
close all

% Sets the paths.
config.path.conn  = '../../data/connectivity/plv_sources_1/';
config.path.str   = '../../data/connectivity/plv_strength_1/';
config.path.patt  = '*.mat';

% Action when the task have already been processed.
config.overwrite  = true;

% Defines the bands to use.
config.bands      = { 'Theta' 'Alpha' 'Low-beta' 'High-beta' 'Beta' 'Gamma' };


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Creates and output folder, if needed.
if ~exist ( config.path.str, 'dir' ), mkdir ( config.path.str ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.conn, config.path.patt ) );

% Goes through each subject.
for file = 1: numel ( files )
    
    % Loads the connectivity data.
    filename   = files ( file ).name;
    [ ~, basename ] = fileparts ( filename );
    
    conndata   = load ( sprintf ( '%s%s', config.path.conn,  basename ), 'subject', 'task', 'stage', 'channel', 'whitener', 'lambda', 'atlas' );
    
    if exist ( sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s.mat', config.path.str, conndata.subject, conndata.task, conndata.stage, conndata.channel, conndata.whitener, conndata.lambda ( 1: end -1 ), conndata.atlas ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'', whitening ''%s'', regularization %s (already calculated).\n', conndata.subject, conndata.task, conndata.channel, conndata.whitener, conndata.lambda );
        continue
    end
    
    fprintf ( 1, 'Working on subject ''%s'', task ''%s'', channel group ''%s'', whitening ''%s'', regularization %s.\n', conndata.subject, conndata.task, conndata.channel, conndata.whitener, conndata.lambda );
    
    fprintf ( 1, '  Loading data.\n' );
    
    conndata   = load ( sprintf ( '%s%s', config.path.conn,  basename ) );
    
    
    fprintf ( 1, '  Loading atlas ''%s''.\n', conndata.atlas );
    
    % Loads the template.
    template   = load ( sprintf ( 'template_%s', conndata.atlas ), 'grid', 'atlas' );
    template   = my_sources2roi ( template );
    
    % Gets the selected areas.
    cortical   = ~ismember ( template.atlas.oarea, [ 0 41 42 71 72 73 74 75 76 77 78 ] );
    
    % Modifies the grid to marke the real inside sources.
    template.grid.inside = ismember ( template.grid.area, find ( cortical ) );
    template.grid.area (:) = 0;
    template.grid.area ( template.grid.inside ) = template.atlas.oarea ( cortical );
    
    
    fprintf ( 1, '  Calculating the strength for each metric and band.\n' );
    
    % Keeps only the requested bands.
    bands  = conndata.band;
    bands  = bands ( ismember ( { bands.name }, config.bands ) );
    
    % Converts the structure into a cell.
    bands  = num2cell ( bands );
    
    % Goes through each band.
    for bindex = 1: numel ( bands )
        
        % Gets the current band.
        band = bands { bindex };
        
        % Selects only the cortical sources.
        band.plv_rms   = band.plv_rms   ( cortical, cortical );
        band.ciplv_rms = band.ciplv_rms ( cortical, cortical );
        
        % Adds the area label.
        band.area      = template.atlas.oarea ( cortical );
        
        % Calculates the strength for each source.
        band.plv_str   = mean ( band.plv_rms, 2, 'omitnan' );
        band.ciplv_str = mean ( band.ciplv_rms, 2, 'omitnan' );
        
        % Removes the unwanted fields.
        band  = rmfield ( band, setdiff ( fieldnames ( band ), { 'name' 'edges' 'length' 'area' 'plv_str' 'ciplv_str' } ) );
        
        % Stores the strength.
        bands { bindex } = band;
    end
    
    % Concatenates all the bands.
    bands = cat ( 1, bands {:} );
    
    
    fprintf ( 1, '  Saving the strength data.\n' );
    
    % Prepares the output.
    strdata           = [];
    strdata.subject   = conndata.subject;
    strdata.task      = conndata.task;
    strdata.stage     = conndata.stage;
    strdata.channel   = conndata.channel;
    strdata.whitener  = conndata.whitener;
    strdata.lambda    = conndata.lambda;
    strdata.atlas     = conndata.atlas;
    strdata.grid      = template.grid;
    strdata.band      = bands;
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s', config.path.str, strdata.subject, strdata.task, strdata.stage, strdata.channel, strdata.whitener, strdata.lambda ( 1: end -1 ), strdata.atlas ), '-struct', 'strdata' );
end
