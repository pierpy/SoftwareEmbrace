function sources = my_mne ( cfg, sourcemodel, data )

% MNE_operator = my_mne ( source_model, noise_cov, source_cov )
% 
% Calculates the inverse MNE operator from the lead field and, if provided,
% the noise and source covariance.

% Based on FieldTrip 20200130 functions:
% * minimumnormestimate by Robert Oostenveld
%
% Based on papers:
% * Dale et al. 2000 Neuron 26.55-67
%   doi: 10.1016/S0896-6273(00)81138-1
% * Liu, Dale & Belliveau 2002 Hum. Brain Mapp. 16.47-62.
%   doi: 10.1002/hbm.10024
% * Lin et al. 2004 Neuroimage 23.582-95.
%   doi: 10.1016/j.neuroimage.2004.04.027
% 
% Based on web-pages:
% * https://mne.tools/0.17/manual/source_localization/inverse.html

% Sets he default parameters.
if ~isfield ( cfg, 'noisecov' ),    cfg.noisecov    = [];    end
if ~isfield ( cfg, 'sourcecov' ),   cfg.sourcecov   = [];    end
if ~isfield ( cfg, 'prewhiten' ),   cfg.prewhiten   = true;  end
if ~isfield ( cfg, 'scalecov' ),    cfg.scalecov    = true;  end
if ~isfield ( cfg, 'lambda' ),      cfg.labmda      = [];    end
if ~isfield ( cfg, 'noiselambda' ), cfg.noiselambda = 0;     end
if ~isfield ( cfg, 'snr' ),         cfg.snr         = [];    end

if ~isfield ( cfg, 'keepfilter' ),  cfg.keepfilter  = true;  end
if ~isfield ( cfg, 'keepmom' ),     cfg.keepmom     = false; end
if ~isfield ( cfg, 'projectmom' ),  cfg.projectmom  = false; end
if ~isfield ( cfg, 'keepcov' ),     cfg.keepcov     = false; end
if ~isfield ( cfg, 'keepnoise' ),   cfg.keepnoise   = false; end

% Gets the options as booleans, if required.
if ischar ( cfg.keepfilter ), cfg.keepfilter = strcmp ( cfg.keepfilter, 'yes' ); end
if ischar ( cfg.keepmom ),    cfg.keepmom    = strcmp ( cfg.keepmom,    'yes' ); end
if ischar ( cfg.projectmom ), cfg.projectmom = strcmp ( cfg.projectmom, 'yes' ); end
if ischar ( cfg.keepcov ),    cfg.keepcov    = strcmp ( cfg.keepcov,    'yes' ); end
if ischar ( cfg.keepnoise ),  cfg.keepnoise  = strcmp ( cfg.keepnoise,  'yes' ); end

if nargin < 3, data = []; end


% Checks the inputs.
if isempty ( cfg.lambda ) && isempty ( cfg.snr )
    error ( 'Either SNR or lambda must be selected.' )
end

% Extracts the lead field.
leadfield   = sourcemodel.leadfield ( sourcemodel.inside );

% Gets the dimensions of the data.
nchan       = numel ( sourcemodel.label );
nsource     = numel ( sourcemodel.inside );
nori        = cellfun ( @(lf) size ( lf, 2 ), leadfield );

% Re-writes the lead field in matrix form.
leadfield   = cat ( 2, leadfield {:} );
leadfield   = double ( leadfield );

% Defines the default covariance matrices.
if ~isempty ( cfg.noisecov )
    noisecov    = double ( cfg.noisecov );
else
    noisecov    = eye ( nchan );
end
if ~isempty ( cfg.sourcecov )
    sourcecov   = double ( cfg.sourcecov );
else
    sourcecov   = speye ( sum ( nori ) );
end


% Generates the whitener, if required.
if cfg.prewhiten
    
    % Gets the PCA decomposition of the covariance matrix.
    [ u, s, ~ ] = svd ( noisecov + cfg.noiselambda * eye ( size ( noisecov ) ) );
    s           = diag ( s );
    sel         = s > s (1) * 1e-12;
%     sel         = s > s (1) * 1e-8;
    
    % Calculates the whitener from the PCA.
    whitener    = diag ( 1 ./ sqrt ( s ( sel ) ) ) * u ( :, sel )';
    whitener    = double ( whitener );
    
    % Pre-whitenes the lead field and the noise covariance.
    leadfield   = whitener * leadfield;
    wnoisecov   = eye ( sum ( sel ) );
    
% Otherwise defines the whitener as the identity matrix.
else
    whitener    = eye ( nchan );
    wnoisecov   = noisecov;
end


% Scales the source covariance matrix, if required.
if cfg.scalecov
    
    % Scales the source covariance so the power is preserved.
    scale       = trace ( leadfield * sourcecov * leadfield' ) / trace ( wnoisecov );
    sourcecov   = sourcecov / scale;
end


% Calculates the lambda parameter form the SNR estimation.
if isempty ( cfg.lambda )
    cfg.lambda = trace ( leadfield * sourcecov * leadfield' ) / ( trace ( wnoisecov ) * snr .^ 2 );
end


if cfg.prewhiten
    
    % Uses the Cholesky factorization to avoid numeric inestabilities.
    sourcecovch = chol ( sourcecov, 'lower' );
    
    % Calculates the regularized inverse using SVD.
    [ u, s, v ] = svd ( leadfield * sourcecovch, 'econ' );
    s           = diag ( s );
    
    % Calculates the minimum norm operator.
    mnefilter   = sourcecovch * ( v * diag ( s ./ ( s .^ 2 + cfg.lambda ) ) * u' );
    
    % Undoes the pre-whitening.
    mnefilter   = mnefilter * whitener;
    
else
    
    % Calculates the minimum norm operator using the explicit formula.
    denom       = leadfield * sourcecov * leadfield' + cfg.lambda .^ 2 * wnoisecov;
%     mnefilter   = ( sourcecov * leadfield' ) / denom;
    mnefilter   = ( sourcecov * leadfield' ) * pinv ( denom );
end


% If required, calculates the moment and the power.
if ~isempty ( data ) && ( cfg.keepmom || cfg.keepcov || cfg.keeppow )
    
    % Gets the number of samples
    nsample     = numel ( data.time );
    
    % Gets the moment from the average time series.
    sourcemom   = mnefilter * data.avg;
    
    % Re-writes the moment as a cell-array.
    sourcemom   = mat2cell ( sourcemom, nori, nsample );
    
    % Calculates the power time series for each source.
    sourcepow   = cellfun ( @(mom) sum ( mom .^ 2, 1 ), sourcemom, 'UniformOutput', false );
    sourcepow   = cat ( 1, sourcepow {:} );
    
    % Calculates the source covariance matrix, if required.
    if cfg.keepcov
        sourcedcov  = cellfun ( @(mom) ( mom * mom' ) / nsample, sourcemom, 'UniformOutput', false );
    end
end


% Re-writes the MNE operator as a cell-array.
mnefilter   = mat2cell ( mnefilter, nori, nchan );

% Calculates the projected noise.
sourcencov  = cellfun ( @(mno) mno * noisecov * mno', mnefilter, 'UniformOutput', false );
sourcenoise = cellfun ( @(n) sum ( diag ( n ) ), sourcencov );


% Initializes the sources structure.
sources              = rmfield ( sourcemodel, setdiff ( fieldnames ( sourcemodel ), { 'pos' 'tri' 'nrm' 'inside' 'unit' 'label' } ) );

% Stores the minimum norm estimator, if required.
if cfg.keepfilter
    sources.filter       = cell ( nsource, 1 );
    sources.filter   ( sources.inside ) = mnefilter;
    sources.filterdimord = '{pos}_ori_chan';
end

% Stores the source-level time series, if required.
if cfg.keepmom
    sources.time         = data.time;
    sources.mom          = cell ( nsource, 1 );
    sources.mom      ( sources.inside ) = sourcemom;
    sources.pow          = nan ( nsource, nsample );
    sources.pow      ( sources.inside, : ) = sourcepow;
end

% Stores the source-level covariance matrices, if required.
if cfg.keepcov
    sources.cov          = cell ( nsource, 1 );
    sources.cov      ( sources.inside ) = sourcedcov;
end

% Stores the source-level noise information, if required.
if cfg.keepnoise
    sources.noise        = nan ( nsource, 1 );
    sources.noise    ( sources.inside ) = sourcenoise;
    sources.noisecov     = cell ( nsource, 1 );
    sources.noisecov ( sources.inside ) = sourcencov;
end
