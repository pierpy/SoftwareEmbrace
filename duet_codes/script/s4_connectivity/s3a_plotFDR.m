clc
clear
close all

% Sets the paths.
config.path.stats = '../../stats/connectivity/plv/';
config.path.figs  = '../../figs/stats/connectivity/fdr/';
config.path.patt  = '*.mat';

config.showthres  = true;


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
    pval1    = statdata.pvalue ( index1, index1 );
    pval2    = statdata.pvalue ( index2, index2 );
    pval12   = statdata.pvalue ( index1, index2 );
    
    
    % Creates the figure.
    figure ( 'Units', 'centimeters', 'Position', [  0.0  1.0 15.0  6.0 ] )
    
    
    % Adds the labels for the participants.
    axes ( 'Units', 'centimeters', 'Position', [  0.0  5.0  5.0  1.0 ] )
    text ( 0, 0, 'Intra-performer 1', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
    xlim ( [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'Position', [  5.0  5.0  5.0  1.0 ] )
    text ( 0, 0, 'Intra-performer 2', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
    xlim ( [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'Position', [ 10.0  5.0  5.0  1.0 ] )
    text ( 0, 0, 'Inter-performers', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
    xlim ( [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    
    % Plots the FDR graph for the first participant.
    axes ( 'Units', 'centimeters', 'Position', [  0.1  0.1  4.8  4.8 ], 'NextPlot', 'add' )
    axis off
    
    % Plots the p-values.
    psort = sort ( pval1 ( isfinite ( pval1 (:) ) ) );
    ax = plot ( psort );
    
    % Plots the random and FDR-corrected random lines.
    plot ( ( 1: numel ( psort ) ) / numel ( psort ), 'k:' )
    plot ( 0.05 * ( 1: numel ( psort ) ) / numel ( psort ), 'r:' )
    
    % Plots the p-threshold, if any.
    if config.showthres && any ( psort < 0.05 * ( 1: numel ( psort ) )' / numel ( psort ) )
        thres = find ( psort < 0.05 * ( 1: numel ( psort ) )' / numel ( psort ), 1, 'last' );
        row = dataTipTextRow ( 'p-value', psort );
        ax.DataTipTemplate.DataTipRows = row;
        datatip ( ax, thres, psort ( thres ), 'FontSize', 7 );
    end
    
    
    % Plots the FDR graph for the second participant.
    axes ( 'Units', 'centimeters', 'Position', [  5.1  0.1  4.8  4.8 ], 'NextPlot', 'add' )
    axis off
    
    % Plots the p-values.
    psort = sort ( pval2 ( isfinite ( pval2 (:) ) ) );
    ax = plot ( psort );
    
    % Plots the random and FDR-corrected random lines.
    plot ( ( 1: numel ( psort ) ) / numel ( psort ), 'k:' )
    plot ( 0.05 * ( 1: numel ( psort ) ) / numel ( psort ), 'r:' )
    
    % Plots the p-threshold, if any.
    if config.showthres && any ( psort < 0.05 * ( 1: numel ( psort ) )' / numel ( psort ) )
        thres = find ( psort < 0.05 * ( 1: numel ( psort ) )' / numel ( psort ), 1, 'last' );
        row = dataTipTextRow ( 'p-value', psort );
        ax.DataTipTemplate.DataTipRows = row;
        datatip ( ax, thres, psort ( thres ), 'FontSize', 7 );
    end
    
    
    % Plots the FDR graph for the inter-participant connectivity.
    axes ( 'Units', 'centimeters', 'Position', [ 10.1  0.1  4.8  4.8 ], 'NextPlot', 'add' )
    axis off
    
    % Plots the p-values.
    psort = sort ( pval12 ( isfinite ( pval12 (:) ) ) );
    ax = plot ( psort );
    
    % Plots the random and FDR-corrected random lines.
    plot ( ( 1: numel ( psort ) ) / numel ( psort ), 'k:' )
    plot ( 0.05 * ( 1: numel ( psort ) ) / numel ( psort ), 'r:' )
    
    % Plots the p-threshold, if any.
    if config.showthres && any ( psort < 0.05 * ( 1: numel ( psort ) )' / numel ( psort ) )
        thres = find ( psort < 0.05 * ( 1: numel ( psort ) )' / numel ( psort ), 1, 'last' );
        row = dataTipTextRow ( 'p-value', psort );
        ax.DataTipTemplate.DataTipRows = row;
        datatip ( ax, thres, psort ( thres ), 'FontSize', 7 );
    end
    
    
    % Saves the figure.
    print ( '-dpng', '-r300', sprintf ( '%s%s_%s_%s_%s.png', config.path.figs, statdata.subject, statdata.comparison, statdata.metric, statdata.bandname ) )
    close
end
