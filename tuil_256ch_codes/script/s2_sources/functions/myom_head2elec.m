function weights = myom_head2elec ( headmodel, elec )


% load ws_danielsson.mat
% elec = sens;



% Gets the electrode positions.
elepos = elec.elecpos;
nele   = size ( elepos, 1 );

% Gets the active surface mesh (the scalp).
mesh   = headmodel.bnd ( end );
npos   = size ( mesh.pos, 1 );
ntri   = size ( mesh.tri, 1 );


% Initializes the weights matrix.
weights = sparse ( nele, npos );

% Goes through each electrode position.
for eindex = 1: nele
    
    % Gets the distance from each element.
    [ dist, weight ] = myom_dpc ( mesh, elec.elecpos ( eindex, : ) );
    
    % Gets the index of closest element.
    [ ~, index ] = min ( dist );
    
    % Gets the nodes of the element.
    nodes = mesh.tri ( index, : );
    
    % Stores the weights.
    weights ( eindex, nodes ) = weight ( index, : );
end

% a = 0;