function source = my_projmne ( source )

% Gets the number of sources and time points.
nsource = numel ( source.inside );
nsample = numel ( source.time );


% Reserves memory for the projection data.
source.phi  = nan ( nsource, 1 );

source.pow  = nan ( nsource, nsample );

% Reserves memory for the z-scored data.
source.zmom = cell ( nsource, 1 );
source.zpow = nan ( nsource, nsample );

% Goes through each inside source.
for sindex = find ( source.inside (:) )'
    
    % Gets the data for the current source alone.
    mom      = source.mom { sindex };
%     noise    = source.noisecov { sindex };
    nrm      = source.nrm ( sindex, : );
    
    % Applies the baseline correction.
    baseline = source.time <= 0;
    mom      = mom - mean ( mom ( :, baseline ), 2 );
    
    % Finds the main direction from the data.
    [ u, ~, ~ ] = svd ( mom * mom' );
    proj     = u ( :, 1 )';
    
    % Projects the data over the main direction.
    mom      = proj * mom;
    phi      = proj * nrm';
%     noise    = proj * noise * proj';
    
    % Stores the projected values.
    source.mom { sindex } = mom;
    source.pow ( sindex, : ) = abs ( mom ) .^ 2;
%     source.noise ( sindex ) = noise;
    source.phi ( sindex ) = phi;
    
    
    % Z-scores the result respect to the baseline.
    zmom     = mom ./ std ( mom ( :, baseline ) );
    zpow     = sum ( zmom .^ 2, 1 );
    
    % Stores the z-scored values.
    source.zmom { sindex } = zmom;
    source.zpow ( sindex, : ) = zpow;
end
