function metadata = my_readmeta ( filename )

% Reads the metadata file.
metadata = readtable ( filename, 'Sheet', 'Metadata', 'PreserveVariableNames', true );

