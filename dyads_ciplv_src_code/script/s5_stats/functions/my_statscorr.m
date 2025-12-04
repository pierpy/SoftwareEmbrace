function pval = my_statscorr ( data, covs, nperms )


% Initializes the input variables.
if nargin < 3, nperms     = 0;  end
if nargin < 2, error ( 'This function requires covariates' ); end

% Gets the number of covariates.
ncov   = size ( covs, 2 );


% If more than one covariate forbids permutations.
if ncov > 1 && nperms > 0
    error ( 'Permutation test can not be performed with more than one covariate.' );
end


% Keeps only the complete data.
comp = all ( isfinite ( covs ), 2 ) & all ( isfinite ( data ), 2 );
data = data ( comp, : );
covs = covs ( comp, : );


% Calculates the original p-value.
pval = stats ( data, covs );


% Checks if permutation statistics are requested.
if nperms
    
    % Initializes the permitation results.
    pperms = zeros ( size ( pval ) );
    
    % Iterates along permutations.
    for pindex = 1: nperms
        
        % Ramdomly re-labels the subjects.
        pcovs = covs ( randperm ( size ( covs, 1 ) ) );
        
        % Calculates the permutation p-value.
        pperm  = stats ( data, pcovs );
        
        % Stores the accumulated p-value.
        pperms = pperms + double ( pperm <= pval );
    end
    
    % Calculates the permutation-corrected p-value.
    pval = pperms / nperms;
end


% Function to get the statistics.
function pval = stats ( data, covs )

% Gets the number of covariates.
ncov   = size ( covs, 2 );

% Gets the N-way ANOVA.
pval   = my_anovan ( data, covs, 'continuous', 1: ncov, 'display', 'off' );

pval   = cat ( 1, pval, pval );
