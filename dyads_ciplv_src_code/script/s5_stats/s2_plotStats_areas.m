clc
clear
close all

% Sets the paths.
config.path.stat = '../../stats/connectivity/plv_areas_window/';
config.path.figs = '../../figs/stats/connectivity/plv_areas_window/';
config.path.patt = '*All*.mat';

% Defines the configuration.
config.qval   = 0.0005;
% config.qval   = 0.10;

% Sets the drawing properties.
config.draw.model           = 'HD';
config.draw.ShowNodesNames  = 'no';
config.draw.transparency    = 0.1;

% Sets the figure width.
% config.width = 8.5;
% config.width = 11.6;
config.width = 17.6;


% Adds the functions folders to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
myft_path


% Creates and output folder, if required.
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ); end


% Lists only the cortical areas.
cortex = setdiff ( 1:90, [ 41 42 71 72 73 74 75 76 77 78 ] );
cortex = cat ( 2, cortex, cortex + 90 );


% Lists the files in the stats folder.
files = dir ( sprintf ( '%s%s', config.path.stat, config.path.patt ) );

% Goes through each file.
for findex = 1: numel ( files )

    % Loads the file.
    statdata = load ( sprintf ( '%s%s', config.path.stat, files ( findex ).name ) );

    fprintf ( 1, 'Working with subject %s, comparison %s, metric %s, band %s (%i-%i Hz).\n', statdata.subject, statdata.comparison, statdata.metric, statdata.bandname, statdata.bandedges );
    
    % Keeps only the cortical areas.
    statdata.label     = statdata.label ( cortex );
    statdata.connmean1 = statdata.connmean1 ( cortex, cortex );
    statdata.connmean2 = statdata.connmean2 ( cortex, cortex );
    statdata.conndiff  = statdata.conndiff ( cortex, cortex );
    statdata.tstat     = statdata.tstat ( cortex, cortex );
    statdata.pvalue    = statdata.pvalue ( cortex, cortex );

    
    % Lists the areas for each participant.
    index1   = ~cellfun ( @isempty, regexp ( statdata.label, '^A_' ) );
    index2   = ~cellfun ( @isempty, regexp ( statdata.label, '^B_' ) );
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
    pthres1  = my_fdr ( pval1, config.qval );
    pthres2  = my_fdr ( pval2, config.qval );
    pthres12 = my_fdr ( pval12, config.qval );





    % Loads the atlas.
    template = load ( 'template_AAL_hyper.mat' );
    atlas    = template.atlas;

    % Keeps only the cortical areas.
    atlas.name  = atlas.name ( cortex );
    atlas.nick  = atlas.nick ( cortex );
    atlas.order = atlas.order ( cortex );
    atlas.pos   = atlas.pos ( cortex, : );
    atlas.conn  = atlas.conn ( cortex, cortex );


    % atlas.name  = atlas.name ( 1: 80 );
    % atlas.nick  = atlas.nick ( 1: 80 );
    % atlas.order = atlas.order ( 1: 80 );
    % atlas.pos   = atlas.pos ( 1: 80, : );
    % atlas.conn  = atlas.conn ( 1: 80, 1: 80 );
    atlas.pos ( index1, 1 ) = atlas.pos ( index1, 1 ) - 0.1;
    atlas.pos ( index2, 1 ) = atlas.pos ( index2, 1 ) + 0.1;





    model = importdata('ModelBrainMNI.mat');
    Faces = model.BrainDavid.Faces;
    Vertices = model.BrainDavid.Vertices / 1e3;

    figure ( 'Units', 'centimeters', 'Position', [  0.00  2.00 10.00  5.00 ] )
    
    hold on
    axis off equal vis3d
    rotate3d

    % Draws the first brain.
    dum = trisurf ( ...
        Faces, ...
        Vertices ( :, 1 ) - 0.1, ...
        Vertices ( :, 2 ), ...
        Vertices ( :, 3 ), ...
        'SpecularStrength', 0.2, ...
        'DiffuseStrength',  0.8, ...
        'AmbientStrength',  0.5, ...
        'FaceLighting',     'phong', ...
        'LineStyle',        'none', ...
        'FaceColor',        [ .75 .75 .75], ...
        'FaceAlpha',        0.1, ...
        'EdgeColor',        'none', ...
        'tag',              'cortex' );

    % Hides the brain in the legend.
    set ( get ( get ( dum, 'Annotation' ), 'LegendInformation' ), 'IconDisplayStyle', 'off' );


    % Draws the second brain.
    dum = trisurf ( ...
        Faces, ...
        Vertices ( :, 1 ) + 0.1, ...
        Vertices ( :, 2 ), ...
        Vertices ( :, 3 ), ...
        'SpecularStrength', 0.2, ...
        'DiffuseStrength',  0.8, ...
        'AmbientStrength',  0.5, ...
        'FaceLighting',     'phong', ...
        'LineStyle',        'none', ...
        'FaceColor',        [ .75 .75 .75], ...
        'FaceAlpha',        0.1, ...
        'EdgeColor',        'none', ...
        'tag',              'cortex' );

    % Hides the brain in the legend.
    set ( get ( get ( dum, 'Annotation' ), 'LegendInformation' ), 'IconDisplayStyle', 'off' );


    % Draws the area centroids.
    scatter3 ( atlas.pos ( :, 1 ), atlas.pos ( :, 2 ), atlas.pos ( :, 3 ), 5, 'fill', 'MarkerFaceColor', [ .4 .4 .4 ] )


    % Draws the links.
    [ i1, i2 ] = find ( pval12 < pthres12 );
    dir = sign ( tstat12 ( pval12 < pthres12 ) );
    i2 = i2 + 80;

    scatter3 ( atlas.pos ( i1, 1 ), atlas.pos ( i1, 2 ), atlas.pos ( i1, 3 ), 25, 'fill', 'MarkerFaceColor', [ .4 .4 .4 ] )
    scatter3 ( atlas.pos ( i2, 1 ), atlas.pos ( i2, 2 ), atlas.pos ( i2, 3 ), 25, 'fill', 'MarkerFaceColor', [ .4 .4 .4 ] )
    x = cat ( 2, atlas.pos ( i1, 1 ), atlas.pos ( i2, 1 ) )';
    y = cat ( 2, atlas.pos ( i1, 2 ), atlas.pos ( i2, 2 ) )';
    z = cat ( 2, atlas.pos ( i1, 3 ), atlas.pos ( i2, 3 ) )';
    plot3 ( x ( :, dir < 0 ), y ( :, dir < 0 ), z ( :, dir < 0 ), 'Color', [ 0.0000 0.4470 0.7410 ] )
    plot3 ( x ( :, dir > 0 ), y ( :, dir > 0 ), z ( :, dir > 0 ), 'Color', [ 0.8500 0.3250 0.0980 ] )
    return


    % Reorders the atlas.
    [ ~, idx ] = sort ( atlas.order );
    order  = zeros ( size ( atlas.order ) );
    order ( idx ) = 1: numel ( atlas.order );


    cfg        = [];
    cfg.width  = config.width;
    cfg.nodes  = atlas.pos ( order, : ) * 1000;
    cfg.groups = statdata.tasks;
    cfg.label  = atlas.nick ( order );

    my_drawNetworkFull ( pval12 ( order, order ) < pthres12, cfg );
    return
    my_drawNetworkFull ( matrix ( order, order ), cfg );
    set ( gcf, 'Name', outname )


    return
    
    % Creates the figure.
    figure ( 'Units', 'centimeters', 'Position', [  0.0  1.0 15.0 16.0 ] )
    
    
    % Adds the labels for the participants.
    axes ( 'Units', 'centimeters', 'Position', [  0.0 15.0  7.5  1.0 ] )
    text ( 0, 0, 'Player 1', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
    xlim ( [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'Position', [  7.5 15.0  7.5  1.0 ] )
    text ( 0, 0, 'Player 2', 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlign', 'middle', 'HorizontalAlign', 'center' )
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
    
    % return
    % Saves the figure.
    print ( '-dpng', '-r300', sprintf ( '%s%s_%s_%s_%s.png', config.path.figs, statdata.subject, statdata.comparison, statdata.metric, statdata.bandname ) )
    close
end
