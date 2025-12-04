function pval = my_stats ( data, group, factors, covars, nperms )


% Initializes the input variables.
if nargin < 5, nperms     = 0;  end
if nargin < 4, covars     = []; end
if nargin < 3, factors    = []; end

% Gets the number of grouping variables, factors, and covariates.
ngroup = size ( group, 2 );
nfact  = size ( factors, 2 ); %#ok<NASGU>
ncov   = size ( covars, 2 ); %#ok<NASGU>


% If more than one grouping variable forbids permutations.
if ngroup > 1 && nperms > 0
    error ( 'Permutation test can not be performed with more than one grouping variable.' );
end


% Calculates the original p-value.
pval = stats ( data, group, factors, covars );


% Checks if permutation statistics are requested.
if nperms
    
    % Initializes the permitation results.
    pperms = zeros ( size ( pval ) );
    
    % Iterates along permutations.
    for pindex = 1: nperms
        
        % Ramdomly re-labels the subjects.
        pgroup = group ( randperm ( numel ( group ) ) );
        
        % Calculates the permutation p-value.
        pperm  = stats ( data, pgroup, factors, covars );
        
        % Stores the accumulated p-value.
        pperms = pperms + double ( pperm <= pval );
    end
    
    % Calculates the permutation-corrected p-value.
    pval = pperms / nperms;
end


% Function to get the statistics.
function pval = stats ( data, group, factors, covars )

% Gets the number of grouping variables and of covariates, if any.
ngroup = size ( group, 2 );
nfact  = size ( factors, 2 );
ncov   = size ( covars, 2 );

% Initializes the input arguments.
inputs    = { 'display', 'off' };

% Rewrites the grouping variables to match anovan inputs.
fullgroup = num2cell ( group, 1 );

% If factors, includes them .
if nfact
    fullgroup = cat ( 2, fullgroup, num2cell ( factors, 1 ) );
end

% If covariates, includes them as continuous grouping variables.
if ncov
    fullgroup = cat ( 2, fullgroup, num2cell ( covars, 1 ) );
    inputs    = cat ( 2, inputs, { 'continuous', ngroup + nfact + ( 1: ncov ) } );
end

% Gets the N-way ANOVA.
[ pval, ~, stats ] = my_anovan ( data, fullgroup, inputs {:} );

% Calculates the pairwise p-value with Tukey's HSD, if needed.
if numel ( fullgroup ) > 1 || numel ( unique ( fullgroup {1} ) ) > 2
    pairwise  = my_multcompare ( stats, 'display', 'off' );
    if size ( pairwise, 1 ) == 1
        pval      = cat ( 1, pval, squeeze ( pairwise ( :, 6, : ) )' );
    else
        pval      = cat ( 1, pval, squeeze ( pairwise ( :, 6, : ) ) );
    end
else
    pval = cat ( 1, pval, pval );
end
