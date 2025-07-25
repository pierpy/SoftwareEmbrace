function source = my_selectOri ( source, threshold )

% If no source covariance does nothing.
if ~isfield ( source, 'cov' )
    return
end


% Extracts the data.
filters       = source.filter;
scovs         = source.cov;


% Goes through each source.
for sindex = 1: numel ( filters )
    
    % Gets the data for the current source.
    filter        = filters { sindex };
    scov          = scovs { sindex };
    
    % Ignores the empty filters.
    if isempty ( filter )
        continue
    end
    
    
    % Calculates the PCA of the orientations.
    [ r, s, ~ ]   = svd ( scov );
    
    
    % Gets the requested cummulative variance.
    if threshold < 1
        cumvar        = cumsum ( diag ( s ) ) / sum ( diag ( s ) );
        ori           = 1: find ( cumvar >= threshold, 1, 'first' );
        r             = r ( ori, : );
        s             = s ( ori, ori );
    end
    
    % Gets the selected number of orientations.
    if threshold >= 1
        ori           = 1: threshold;
        r             = r ( ori, : );
        s             = s ( ori, ori );
    end
    
    
    % Rotates the filter according to the PCA.
    filter        = r * filter;
    scov          = s;
    
    % Stores the filter.
    filters { sindex } = filter;
    scovs   { sindex } = scov;
end

% Stores the data.
source.filter = filters;
source.cov    = scovs;
