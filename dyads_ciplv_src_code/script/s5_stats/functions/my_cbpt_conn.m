function stat = my_cbpt_conn ( data, group, covariates, anova, parametric, nperm )

% Initializes the input variables.
if nargin < 6, nperm      = 0;  end
if nargin < 5, parametric = 1;  end
if nargin < 4, anova      = 0;  end
if nargin < 3, covariates = []; end

config.calpha = 0.01;


% Gets the number of grouping variables and of covariates, if any.
ngroup = size ( group, 2 );
ncov   = size ( covariates, 2 );
nnod   = size ( data, 2 );

% If there are covariates the test must be ANOVA.
if ncov && ~parametric
    error ( 'Only parametric ANOVA can include covariates.' );
end

if ncov && ~anova
    warning ( 'Forcing the test to ANOVA to include covariates.' );
    anova = 1;
end

% If more than one grouping variable allows only ANOVA.
if ngroup > 1 && ~parametric
    error ( 'Only parametric ANOVA can use more than one grouping variable.' );
end

if ngroup > 1 && ~anova
    warning ( 'Forcing the test to ANOVA to allow more than one grouping variable.' );
    anova = 1;
end

% If more than one grouping variable forbids permutations.
if ngroup > 1
    error ( 'Permutation test can not be performed with more than one grouping variable.' );
end


% Removes the repeated values in simmetrical measures.
[ data, ~, dindex ] = unique ( data ( :, : )', 'rows' );
data = data';


% Calculates the statistics.
[ ~, pval, ~, stats ] = ttest2 ( data ( group == 2, : ), data ( group == 1, : ) );

% Restores the original shape of the data.
pval = pval ( dindex );
tval = stats.tstat ( dindex );

pmat = reshape ( pval, nnod, nnod );
tmat = reshape ( tval, nnod, nnod );

pmat ( eye ( nnod ) == 1 ) = 1;
tmat ( eye ( nnod ) == 1 ) = 0;


% Extracts the clusters.
matrix = tmat .* ( pmat <= config.calpha );

clus   = my_connclusters ( matrix );
nclus  = size ( clus, 2 );

% Calculates the statistic for each cluster.
tclus  = zeros ( nclus, 1 );
for cindex = 1: nclus
    tclus ( cindex ) = nansum ( nansum ( matrix ( clus ( :, cindex ) ~= 0, clus ( :, cindex ) ~= 0 ) ) ) / 2;
end

% Sorts the clusters by t-value.
[ ~, index ] = sort ( abs ( tclus ), 'descend' );
clus   = clus ( :, index );
tclus  = tclus ( index );

% Labels the connectivity matrix.
mat = double ( matrix ~= 0 );
for cindex = 1: nclus
    cmember = clus ( :, cindex ) ~= 0;
    mat ( cmember, : ) = cindex * ( matrix ( cmember, : ) ~= 0 );
    mat ( :, cmember ) = cindex * ( matrix ( :, cmember ) ~= 0 );
end


% Initializes the vector of cluster p-values.
pclus = zeros ( nclus, 1 );

% Calculates the randome permutations' clusters.
for iindex = 1: nperm
    
    % Permutes the labels.
    giter = group ( randperm ( size ( group, 1 ) ) );
    
    % Calculates the statistics.
    [ ~, pval, ~, stats ] = ttest2 ( data ( giter == 1, : ), data ( giter == 2, : ) );
    
    % Restores the original shape of the data.
    pval = pval ( dindex );
    tval = stats.tstat ( dindex );
    
    pperm = reshape ( pval, nnod, nnod );
    tperm = reshape ( tval, nnod, nnod );
    
    pperm ( eye ( nnod ) == 1 ) = 1;
    tperm ( eye ( nnod ) == 1 ) = 0;
    
    
    % Extracts the clusters.
    matrix = tperm .* ( pperm <= config.calpha );
    
    citer  = my_connclusters ( matrix );
    nciter = size ( citer, 2 );
    
    % If no cluster, ignores the permutation.
    if nciter == 0
        continue
    end
    
    % Calculates the statistic for each cluster.
    pciter = zeros ( nciter, 1 );
    for cindex = 1: nciter
        pciter ( cindex ) = nansum ( nansum ( matrix ( citer ( :, cindex ) ~= 0, citer ( :, cindex ) ~= 0 ) ) ) / 2;
    end
    
    % Compares the bigger permutation cluster with the original clusters.
    thres = max ( abs ( pciter ) );
    pclus = pclus + double ( abs ( tclus ) < thres );
end

% Divides by the number of permutations.
pclus = pclus / ( nperm + 1 );

% % Corrects the two tails.
% pclus = min ( 2 * pclus, 1 );


% Sets the output.
stat.pval  = pval;
stat.mat   = mat;
stat.clus  = clus;
stat.tclus = tclus;
stat.pclus = pclus;
