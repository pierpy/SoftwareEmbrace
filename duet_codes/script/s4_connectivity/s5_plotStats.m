clc
clear
close all

% Sets the paths.
config.path.stats = '../../stats/connectivity/plv_sessions/';
config.path.figs  = '../../figs/stats/connectivity/plv_sessions/';
config.path.patt  = '*_ciplv_*.mat';


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Creates and output folder, if required.
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ); end


% Loads the layout.
sensmeta = load ( 'ant64dry_hyper.mat' );


% Splits the layout.
layout   = sensmeta.layout;
index1   = ~cellfun ( @isempty, regexp ( sensmeta.label, '_1$' ) );
index2   = ~cellfun ( @isempty, regexp ( sensmeta.label, '_2$' ) );

% Displaces the channel positions.
layout.pos ( index1, 1 ) = layout.pos ( index1, 1 ) - 0.6;
layout.pos ( index2, 1 ) = layout.pos ( index2, 1 ) + 0.6;

% Duplicates and displaces the outline of the head.
outline  = layout.outline;
outline1 = cellfun ( @(pos) pos - [ 0.6 0.0 ], outline, 'UniformOutput', false );
outline2 = cellfun ( @(pos) pos + [ 0.6 0.0 ], outline, 'UniformOutput', false );
outline  = cat ( 2, outline1, outline2 );
layout.outline = outline;


% Lists the files in the stats folder.
files = dir ( sprintf ( '%s%s', config.path.stats, config.path.patt ) );

% Goes through each file.
for findex = 1: numel ( files )
    
    % Loads the file.
    statdata = load ( sprintf ( '%s%s', config.path.stats, files ( findex ).name ) );
    
    fprintf ( 1, 'Working with subject %s, comparison %s, metric %s, band %s (%i-%i Hz).\n', statdata.subject, statdata.comparison, statdata.metric, statdata.bandname, statdata.bandedges );
    
    
    % Lists the channels for each participant.
    index1   = ~cellfun ( @isempty, regexp ( statdata.label, '_1$' ) );
    index2   = ~cellfun ( @isempty, regexp ( statdata.label, '_2$' ) );
    label1   = statdata.label ( index1 );
    label2   = statdata.label ( index2 );
    
    % Gets the inter- and inta-participant matrices.
    tstat1   = statdata.tstat ( index1, index1 );
    pval1    = statdata.pvalue ( index1, index1 );
    tstat2   = statdata.tstat ( index2, index2 );
    pval2    = statdata.pvalue ( index2, index2 );
    tstat12  = statdata.tstat ( index1, index2 );
    pval12   = statdata.pvalue ( index1, index2 );
    
    % Determines the FDR-corrected threshold for the p-values.
    pthres1  = my_fdr ( pval1 );
    pthres2  = my_fdr ( pval2 );
    pthres12 = my_fdr ( pval12 );
    
    
    % Creates the figure.
    figure ( 'Units', 'centimeters', 'Position', [  0.0  1.0 15.0 16.0 ] )
    
    
    % Adds the labels for the participants.
    axes ( 'Units', 'centimeters', 'Position', [  0.0 15.0  7.5  1.0 ] )
    text ( 0, 0, 'Performer 1', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
    xlim ( [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'Position', [  7.5 15.0  7.5  1.0 ] )
    text ( 0, 0, 'Performer 2', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
    xlim ( [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    
    % Plots the base layout for the intra-participant connectivity.
    axes ( 'Units', 'centimeters', 'Position', [  0.1  7.6 14.8  7.3 ] )
    ft_plot_layout ( layout, ...
        'chanindx', find ( index1 | index2 ), ...
        'mask', false, ...
        'point', true, ...
        'box', false, ...
        'label', false )
    hold on
    
    
    % Gets the intra-participant significant links for participant 1.
    hits      = find ( pval1 <= pthres1 );
    [ hit1, hit2 ] = ind2sub ( size ( pval1 ), hits );
    dir       = sign ( tstat1 ( hits ) );
    
    % Gets the position of the electrodes in the layout.
    hitind1   = my_matchstr ( layout.label, label1 ( hit1 ) );
    hitpos1   = layout.pos ( hitind1, : );
    hitind2   = my_matchstr ( layout.label, label1 ( hit2 ) );
    hitpos2   = layout.pos ( hitind2, : );
    
    % Plots the significant links.
    hitposx   = cat ( 1, hitpos1 ( :, 1 )', hitpos2 ( :, 1 )' );
    hitposy   = cat ( 1, hitpos1 ( :, 2 )', hitpos2 ( :, 2 )' );
    
    plot ( hitposx ( :, dir < 0 ), hitposy ( :, dir < 0 ), 'Color', '#0072BD' )
    plot ( hitposx ( :, dir > 0 ), hitposy ( :, dir > 0 ), 'Color', '#D95319' )
    
    
    % Gets the intra-participant significant links for participant 2.
    hits      = find ( pval2 <= pthres2 );
    [ hit1, hit2 ] = ind2sub ( size ( pval2 ), hits );
    dir       = sign ( tstat2 ( hits ) );
    
    % Gets the position of the electrodes in the layout.
    hitind1   = my_matchstr ( layout.label, label2 ( hit1 ) );
    hitpos1   = layout.pos ( hitind1, : );
    hitind2   = my_matchstr ( layout.label, label2 ( hit2 ) );
    hitpos2   = layout.pos ( hitind2, : );
    
    % Plots the significant links.
    hitposx   = cat ( 1, hitpos1 ( :, 1 )', hitpos2 ( :, 1 )' );
    hitposy   = cat ( 1, hitpos1 ( :, 2 )', hitpos2 ( :, 2 )' );
    
    plot ( hitposx ( :, dir < 0 ), hitposy ( :, dir < 0 ), 'Color', '#0072BD' )
    plot ( hitposx ( :, dir > 0 ), hitposy ( :, dir > 0 ), 'Color', '#D95319' )
    
    
    % Plots the base layout for the inter-participant connectivity.
    axes ( 'Units', 'centimeters', 'Position', [  0.1  0.1 14.8  7.3 ] )
    ft_plot_layout ( layout, ...
        'chanindx', find ( index1 | index2 ), ...
        'mask', false, ...
        'point', true, ...
        'box', false, ...
        'label', false )
    hold on
    
    % Gets the inter-participant significant links.
    hits      = find ( pval12 <= pthres12 );
    [ hit1, hit2 ] = ind2sub ( size ( pval12 ), hits );
    dir       = sign ( tstat12 ( hits ) );
    
    % Gets the position of the electrodes in the layout.
    hitind1   = my_matchstr ( layout.label, label1 ( hit1 ) );
    hitpos1   = layout.pos ( hitind1, : );
    hitind2   = my_matchstr ( layout.label, label2 ( hit2 ) );
    hitpos2   = layout.pos ( hitind2, : );
    
    % Plots the significant links.
    hitposx   = cat ( 1, hitpos1 ( :, 1 )', hitpos2 ( :, 1 )' );
    hitposy   = cat ( 1, hitpos1 ( :, 2 )', hitpos2 ( :, 2 )' );
    
    plot ( hitposx ( :, dir < 0 ), hitposy ( :, dir < 0 ), 'Color', '#0072BD' )
    plot ( hitposx ( :, dir > 0 ), hitposy ( :, dir > 0 ), 'Color', '#D95319' )
    
    
    % Saves the figure.
    print ( '-dpng', '-r300', sprintf ( '%s%s_%s_%s_%s.png', config.path.figs, statdata.subject, statdata.comparison, statdata.metric, statdata.bandname ) )
    close
end
