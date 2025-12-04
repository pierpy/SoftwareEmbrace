function clusters = my_connclusters ( matrix )


% Gets the size of the data.
nnod       = size ( matrix, 1 );

% Makes sure that the matrix is binary.
if ~islogical ( matrix )
    posmatrix  = matrix > 0;
    negmatrix  = matrix < 0;
else
    posmatrix  = matrix;
    negmatrix  = zeros ( nnod );
end

% Makes sure that the matrix is simmetric and the diagonal is active.
posmatrix  = double ( posmatrix | posmatrix | eye ( nnod ) > 0 );
negmatrix  = double ( negmatrix | negmatrix | eye ( nnod ) > 0 );

% Closes the clusters.
posmatrix  = posmatrix ^ ( nnod / 2 ) > 0;
negmatrix  = negmatrix ^ ( nnod / 2 ) > 0;

% Removes the clusters with only one member.
posmatrix  = posmatrix ( sum ( posmatrix, 2 ) > 1, : );
negmatrix  = negmatrix ( sum ( negmatrix, 2 ) > 1, : );

% Gets the list of clusters.
poscluster = unique ( posmatrix, 'rows' )';
negcluster = unique ( negmatrix, 'rows' )';

% Sorts the clusters.
[ ~, ord ] = sort ( sum ( poscluster, 1 ), 'descend' );
poscluster = poscluster ( :, ord );
[ ~, ord ] = sort ( sum ( negcluster, 1 ), 'descend' );
negcluster = negcluster ( :, ord );

% Joins the positive and negative clusters.
clusters  = cat ( 2, poscluster, -negcluster );
