function pthres = my_fdr ( pval, qval )

if nargin < 2
    qval   = 0.05;
end


% Sorts the p-values.
psort  = sort ( pval (:) );

% Generates the FDR expected list of sorted p-values.
prand  = qval * ( 1: numel ( psort ) )' / numel ( psort );

% Compares the sorted p-values with the FDR expected graph.
thres  = find ( psort <= prand, 1, 'last' );

% Determines the threshold p-value.
pthres = psort ( thres );

% If no hits, sets the threshold to 0.
pthres = max ( cat ( 1, 0, pthres ), [], 1 );
