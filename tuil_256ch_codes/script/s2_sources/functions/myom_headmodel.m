function headmodel = myom_headmodel ( headmodel, expflag )
% Based on FiedTrip functions:
% * ft_headmodel_openmeeg

global my_silent
my_silent = ~isempty ( my_silent ) && my_silent;


% Adds OpenMEEG to the path.
ft_hastoolbox ( 'openmeeg', 1, 1 );

% Checks the integrity of the OpenMEEG binaries.
om_checkombin;

% Checks the geometrical definition of the head model.
headmodel = myom_check_headmodel ( headmodel );

% Checks the input.
if nargin < 2, expflag = false; end


% Sets the temporal base name.
basename = tempname;

% Writes the geometry files.
myom_write_geometry ( basename, headmodel )


% Calculates the head matrix.
if ~my_silent
    status = system ( sprintf ( 'om_assemble -hm "%s.geom" "%s.cond" "%s_hm.mat"\n', basename, basename, basename ) );
else
    [ status, output ] = system ( sprintf ( 'om_assemble -hm "%s.geom" "%s.cond" "%s_hm.mat"\n', basename, basename, basename ) );
end

% Checks for the completion of the execution.
if status ~= 0
    if my_silent, fprintf ( 1, '%s', output ); end
    fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );
    
    % Removes all the temporal files and exits.
    delete ( sprintf ( '%s*', basename ) );
    return
end

% Recovers the calculated head model matrix.
headmodel.hm = importdata ( sprintf ( '%s_hm.mat', basename ) );

% Expands the matrix, if requested.
if expflag
    headmodel.hm = myom_struct2sym ( headmodel.hm );
end

% Removes all the temporal files.
delete ( sprintf ( '%s*', basename ) );


function matrix = myom_struct2sym ( structure )

% Generates a matrix of the given size.
matrix = zeros ( structure.size );

% Copies the upper diagonal from the structure data to the matrix.
matrix ( triu ( true ( structure.size ) ) ) = structure.data;

% Fills the lower diagonal transposing and removing the diagonal.
matrix = matrix + matrix' - diag ( diag ( matrix ) );
