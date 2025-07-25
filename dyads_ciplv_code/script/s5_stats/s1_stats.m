clc
clear
% close all

% Defines the paths to the data.
config.path.conn = '../../data/connectivity/plv_window/';
config.path.stat = '../../stats/connectivity/plv_window/';
config.path.patt = '*.mat';

% Defines the configuration.
config.band   = 'Alpha';
config.metric = 'ciplv';
config.qval   = 0.1;

% Defines the conditions to compare.
config.cond (1).label = 'CoopEasy';
config.cond (1).hits  = [ NaN 1 1 NaN; NaN 1 2 NaN ];
config.cond (1).hits  = [ NaN 1 1 NaN ];
config.cond (2).label = 'CompEasy';
config.cond (2).hits  = [ NaN 2 1 NaN; NaN 2 2 NaN ];
config.cond (2).hits  = [ NaN 2 1 NaN ];
config.cond (1).label = 'CoopAll';
config.cond (1).hits  = [ NaN 1 1 NaN; NaN 1 2 NaN ];
config.cond (2).label = 'CompAll';
config.cond (2).hits  = [ NaN 2 1 NaN; NaN 2 2 NaN ];
config.cond (1).label = 'CoopHard';
config.cond (1).hits  = [ NaN 1 2 NaN ];
config.cond (2).label = 'CompHard';
config.cond (2).hits  = [ NaN 2 2 NaN ];


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions', pwd ) );

% Adds, if required, FieldTrip to the path.
myft_path


% Creates the output folder, if required.
if ~exist ( config.path.stat, 'dir' ), mkdir ( config.path.stat ); end


% Lists the files.
files = dir ( sprintf ( '%s%s', config.path.conn, config.path.patt ) );

% Reserves memory for the data and metadata.
conndatas = cell ( numel ( files ), 1 );
rawconns  = cell ( numel ( files ), 1 );
metadatas = cell ( numel ( files ), 1 );

% Goes through each file.
for findex = 1: numel ( files )
    
    % Loads the file.
    conndata = load ( sprintf ( '%s%s', config.path.conn, files ( findex ).name ) );
    
    % Gets the metadata for this file.
    metadata = conndata.trialinfo.trialdef;
    
    % Adds the file to the metadata.
    metadata ( :, end + 1 ) = findex; %#ok<SAGROW>
    
    % Gets the desired band.
    bindex   = strcmp ( { conndata.band.name }, config.band );
    banddata = conndata.band ( bindex );

    % Gets the connectivity matrix.
    rawconn  = banddata.( config.metric );
    
    % Flattens the connectivity matrix.
    rawconn  = permute ( rawconn, [ 3 1 2 ] );
    % rawconn  = reshape ( rawconn, size ( rawconn, 1 ), [] );

    % Replaces the diagonal by 0s.
    diags    = eye ( numel ( banddata.label ) ) == 1;
    rawconn ( :, diags ) = 0;
    
    % % Keeps only the lower triangular.
    % hits     = tril ( true ( 128 ), -1 );
    % hits     = ( 1: 128 )' > 64 & ( 1: 128 ) < 65;
    % rawconn  = rawconn ( :, hits (:) );
    
    % Stores the data and metadata for this file.
    metadatas { findex } = metadata;
    rawconns  { findex } = rawconn;
    conndatas { findex } = rmfield ( conndata, 'band' );
end

% Concatenates all the data and metadata.
metadatas = cat ( 1, metadatas {:} );
rawconns  = cat ( 1, rawconns  {:} );
conndatas = cat ( 1, conndatas {:} );


% Gets the valid combinations of identifiers.
valid     = cat ( 1, config.cond.hits );
valid     = permute ( valid, [ 3 2 1 ] );

% Identifies the valid conditions.
hits      = metadatas == valid;
hits      = hits | isnan ( valid );
hits      = any ( all ( hits, 2 ), 3 );

% Keeps only the valid conditions.
metadatas = metadatas ( hits, : );
rawconns  = rawconns  ( hits, :, : );


% Builds the design matrix.
design    = nan ( size ( metadatas, 1 ), 1 );

% Goes through each condition.
for cindex = 1: numel ( config.cond )
    
    % Gets the combination of identifiers for this condition.
    valid     = config.cond ( cindex ).hits;
    valid     = permute ( valid, [ 3 2 1 ] );

    % Identifies the valid conditions.
    hits      = metadatas == valid;
    hits      = hits | isnan ( valid );
    hits      = any ( all ( hits, 2 ), 3 );
    
    % Adds the condition to the design matrix.
    design ( hits, 1 ) = cindex;
end


% Performs the comparison.
[ pvalue, ~, stats, ~ ] = my_anovan ( rawconns, design, 'display', 'off' );

% Gets the pairwise statistics.
[ comparison, means, ~, ~, tpair ] = my_multcompare ( stats, 'display', 'off' );

% Gets the p-values in the right shape.
ppair = comparison ( :, 6, : );
ppair = permute ( ppair, [ 1 3 2 ] );

% Gets the means in the right shape.
means = means ( :, 1, : );
means = permute ( means, [ 1 3 2 ] );

% Rewrites all the pair statistics in matrix form.
tpair = reshape ( tpair, [], numel ( banddata.label ), numel ( banddata.label ) );
ppair = reshape ( ppair, [], numel ( banddata.label ), numel ( banddata.label ) );
means = reshape ( means, [], numel ( banddata.label ), numel ( banddata.label ) );
tpair = permute ( tpair, [ 2 3 1 ] );
ppair = permute ( ppair, [ 2 3 1 ] );
means = permute ( means, [ 2 3 1 ] );


% Prepares the output.
statdata  = [];
statdata.subject    = strjoin ( { conndatas.subject }, '+' );
statdata.comparison = strjoin ( { config.cond.label }, '-v-' );
statdata.tasks      = { config.cond.label };
statdata.bandname   = banddata.name;
statdata.bandedges  = banddata.edges;
statdata.metric     = config.metric;
statdata.label      = banddata.label;
statdata.connmean1  = means ( :, :, 1 );
statdata.connmean2  = means ( :, :, 2 );
statdata.conndiff   = means ( :, :, 1 ) - means ( :, :, 2 );
statdata.tstat      = tpair;
statdata.pvalue     = ppair;


save ( '-v6', sprintf ( '%s%s_%s_%s_%s.mat', config.path.stat, statdata.subject, statdata.comparison, statdata.metric, statdata.bandname ), '-struct', 'statdata' )

% psort     = sort ( pvalue );
% fdrthr    = ( 1: numel ( pvalue ) ) / numel ( pvalue ) * config.qval;
% pthr      = psort ( find ( psort <= fdrthr, 1, 'last' ) );
% if isempty ( pthr )
%     pthr = 0;
% end
% 
% 
% plot ( psort )
% hold on
% plot ( fdrthr, ':r' )
% title ( sprintf ( 'p-threshold: %.4f (%i)', pthr, sum ( pvalue <= pthr ) ) )
