clc
clear
close all


% Loads the image.
img     = imread ( 'Layout Eego.png' );

% Removes the reference and the ground.
img     = img ( :, 201: end, : );

% Flips the Y dimension.
img     = flipud ( img );

% Binarizes the image.
imgbw   = img ( :, :, 1 ) < 250;

% Erodes and dilates the image to separate the electrode positions.
imgbw   = imerode ( imgbw, strel ( 'disk', 8 ) );
imgbw   = imopen ( imgbw, strel ( 'disk', 8 ) );

% Identifies the electrodes.
imgelec = bwlabel ( imgbw );
items   = unique ( imgelec ( imgelec (:) ~= 0 ) );


% Reserves memory for the electrode positions.
elecpos = zeros ( numel ( items ), 1 );

% Goes through each electrode.
for eindex = 1: size ( elecpos, 1 )
    
    % Finds the electrode.
    [ y, x ] = find ( imgelec == items ( eindex ) );
    
    % Gets its center of masses.
    elecpos ( eindex, 1 ) = mean ( x );
    elecpos ( eindex, 2 ) = mean ( y );
end

% % Plots the result.
% figure
% imshow ( img )
% hold on
% set ( gca, 'YDir', 'normal' )
% scatter ( elecpos ( :, 1 ), elecpos ( :, 2 ), 'filled' )
% axis equal


% % Initializes the electrode labels.
% labels  = cell ( numel ( items ), 1 );
% 
% % Goes through each electrode.
% for eindex = 1: numel ( labels )
%     
%     % Marks the current electrode in red.
%     scatter ( elecpos ( eindex, 1 ), elecpos ( eindex, 2 ), 'red', 'filled', 'Tag', 'active' )
%     
%     % Asks for the label.
%     label   = mydlg_inputdlg ( 'What is the label of the current electrode?', 'Label' );
%     
%     % If empty, cancels.
%     if isempty ( label ), break, end
%     
%     % Deletes the electrode.
%     delete ( findall ( gca, 'Tag', 'active' ) )
%     
%     % Stores the label.
%     labels ( eindex ) = label;
% end

% Sets the previously defined labels.
labels  = { ...
    '3LD', '2LD', '3LC', '4LC', '2LC', '1LD', '3LB', '4LD', '4LB', '2LB', '5LC', '1LC', '5LB', '2LA', '1LA', '3LA', ...
    '1LB', '10L', '9L', '8L', '7L', '6L', '5L', '1L', '4L', '3L', '2L', '6Z', '8Z', '5Z', '9Z', '7Z', ...
    '1Z', '0Z', '4Z', '2Z', '3Z', '2R', '6R', '3R', '5R', '1R', '7R', '4R', '8R', '9R', '10R', '1RB', ...
    '3RA', '1RA', '2RA', '5RB', '1RC', '2RB', '5RC', '4RB', '3RB', '4RD', '1RD', '2RC', '4RC', '3RC', '3RD', '2RD' };


% % Plots the result.
% figure
% imshow ( img )
% hold on
% set ( gca, 'YDir', 'normal' )
% text ( elecpos ( :, 1 ), elecpos ( :, 2 ), labels, 'VerticalAlign', 'top', 'HorizontalAlign', 'center', 'BackgroundColor', [ 1.0 1.0 1.0 0.8 ], 'FontSize', 14, 'FontWeight', 'bold' )
% scatter ( elecpos ( :, 1 ), elecpos ( :, 2 ), 'filled' )
% axis equal

%%
% Gets the horizontal center of the layout definition.
hcenter = mean ( elecpos ( :, 1 ) );

% Finds the Z electrodes.
hits    = regexp ( labels, '^[\d]+Z$' );
zeroes  = ~cellfun ( @isempty, hits );

% Sets the X coordinate to the center.
elecpos ( zeroes, 1 ) = hcenter;



% Gets the lateralized electrodes.
hits    = regexp ( labels, '^[\d]+L[ABCD]*$' );
lefts   = ~cellfun ( @isempty, hits );
llabel  = labels ( lefts );
rlabel  = strrep ( llabel, 'L', 'R' );

% Goes through each laterlaized electrode.
for eindex = 1: numel ( llabel )
    
    % Gets the positions.
    lpos    = elecpos ( strcmp ( labels, llabel ( eindex ) ), : );
    rpos    = elecpos ( strcmp ( labels, rlabel ( eindex ) ), : );
    
    % The X coordinate is the average distance to the center.
    ldist   = abs ( lpos (1) - hcenter );
    rdist   = abs ( rpos (1) - hcenter );
    xdist   = ( ldist + rdist ) / 2;
    
    % The Y coordinate is the average of both coordinates.
    ypos    = ( lpos (2) + rpos (2) ) / 2;
    
    % Sets the simetrized positions.
    lpos    = [ ( hcenter - xdist ) ypos ];
    rpos    = [ ( hcenter + xdist ) ypos ];
    
    % Stores the positions.
    elecpos ( strcmp ( labels, llabel ( eindex ) ), : ) = lpos;
    elecpos ( strcmp ( labels, rlabel ( eindex ) ), : ) = rpos;
end

% Plots the result.
figure
imshow ( img )
hold on
set ( gca, 'YDir', 'normal' )
text ( elecpos ( :, 1 ), elecpos ( :, 2 ), labels, 'VerticalAlign', 'top', 'HorizontalAlign', 'center', 'BackgroundColor', [ 1.0 1.0 1.0 0.8 ], 'FontSize', 14, 'FontWeight', 'bold' )
scatter ( elecpos ( :, 1 ), elecpos ( :, 2 ), 'filled' )
axis equal


% Creates a table with the layout positions.
layout       = cell2table ( labels (:), 'VariableNames', { 'Label' } );
layout.X     = elecpos ( :, 1 );
layout.Y     = elecpos ( :, 2 );

% Saves the table.
writetable ( layout, 'ant64dry_layout.xlsx' );
