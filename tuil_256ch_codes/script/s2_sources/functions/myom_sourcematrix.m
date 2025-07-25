function headmodel = myom_sourcematrix ( headmodel, grid )

% Based on FiedTrip functions:
% * ft_leadfield_openmeeg by Daniel D.E. Wong, Sarang S. Dalal
%
% Based on the OpenMEEG functions:
% * openmeeg_dsm by Alexandre Gramfort
% * openmeeg_megm by Emmanuel Olivi

global my_silent
my_silent  = ~isempty ( my_silent ) && my_silent;


% Adds OpenMEEG to the path.
ft_hastoolbox ( 'openmeeg', 1, 1 );

% Checks the integrity of the OpenMEEG binaries.
om_checkombin;

% Checks the geometrical definition of the head model.
headmodel = myom_check_headmodel ( headmodel );

% Sanitizes the grid.
grid      = my_fixgrid ( grid );


% Sets the temporal base name.
basename  = tempname;

% Writes the dipole file.
myom_write_dipoles ( basename, grid )

% Writes the geometry files.
myom_write_geometry ( basename, headmodel )

% Determines the domain where the sources are located.
srctissue = headmodel.tissue { headmodel.source };


% Calculates the dipoles matrix.
if ~my_silent
    status = system ( sprintf ( 'om_assemble -dsm "%s.geom" "%s.cond" "%s_dip.bin" "%s_dsm.mat" "%s"\n', basename, basename, basename, basename, srctissue ) );
else
    [ status, output ] = system ( sprintf ( 'om_assemble -dsm "%s.geom" "%s.cond" "%s_dip.bin" "%s_dsm.mat" "%s"\n', basename, basename, basename, basename, srctissue ) );
end

% Checks for the completion of the execution.
if status ~= 0
    if my_silent, fprintf ( 1, '%s', output ); end
    fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );

    % Removes all the temporal files and exits.
    delete ( sprintf ( '%s*', basename ) );
    return
end

% Recovers the calculated dipoles model matrix.
headmodel.dsm = importdata ( sprintf ( '%s_dsm.mat', basename ) );


% If the headmodel matrix is present calculates hm_dsm.
if isfield ( headmodel, 'hm' )
    
    if ~my_silent, fprintf ( 1, 'Calculating inv ( hm ) * dsm.\n' ); end
    
    % Transforms the head model matrix to a symmetric matrix.
    if isstruct ( headmodel.hm )
        headmodel.hm = myom_struct2sym ( headmodel.hm );
    end
    
    % Calculates inv ( hm ) * dsm.
    headmodel.hm_dsm = headmodel.hm \ headmodel.dsm;
    
    % Deletes the head model and sources matrices to save memory.
    headmodel = rmfield ( headmodel, { 'hm' 'dsm' } );
    
elseif isfield ( headmodel, 'ihm' )
    
    if ~my_silent, fprintf ( 1, 'Calculating inv ( hm ) * dsm.\n' ); end
    
    % Calculates inv ( hm ) * dsm.
    headmodel.hm_dsm = headmodel.ihm * headmodel.dsm;
    
    % Deletes the sources matrix to save memory.
    headmodel = rmfield ( headmodel, 'dsm' );
end

% Removes all the temporal files.
delete ( sprintf ( '%s*', basename ) );


% Stores the grid with the head model.
headmodel.grid = grid;



function matrix = myom_struct2sym ( structure )

% Generates a matrix of the given size.
matrix = zeros ( structure.size );

% Copies the upper diagonal from the structure data to the matrix.
matrix ( triu ( true ( structure.size ) ) ) = structure.data;

% Fills the lower diagonal transposing and removing the diagonal.
matrix = matrix + matrix' - diag ( diag ( matrix ) );
