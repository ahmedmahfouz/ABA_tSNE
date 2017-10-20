%% Make the plots for the mouse tsne coronal view (PCA sweep)
clear variables
close all

% Load additional color data
ids = {'9861','10021','12876','14380','15496','15697'};
markersize = 30;
figsubnums = 'abcdef';

%% Figure a-f
% t-SNE mapping colored by region with 10 components
% for idnr = 1:length(ids)
idnr = 1;
% Load data
load(['results_human_average_pca_sweep/MappedHumanGenesAverage2D' ids{idnr} 'id300'],'mappedX','color_RGB');
datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' ids{idnr} '/'];
load([datadir 'sample']);

% Plot figure
nsamples = size(mappedX,1);
figure(1)
clf
hold on
for snr = 1:nsamples
    %     if ~isempty(colors.voxel.color_HEX{snr})
    %         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
    %     else
    %         colorRGB = [0 0 0];
    %     end
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',color_RGB(snr,:),'Markersize',markersize);
    %     text(mappedX(snr,1),mappedX(snr,2),sample.structure_name{snr});
end

% Format the plot
axis off equal tight
set(gcf,'Color','white')
dpi = 300;
figsizeinch = [7 7]/2.54;
set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
% saveas(1,['Figure5' figsubnums(idnr) '.png']);
% end

%% Figure g-l
% t-SNE mapping colored by region with 10 components
figsubnums2 = 'ghijkl';
for idnr = 1:length(ids)
    % idnr = 1;
    % Load data
    load(['results_human_average_pca_sweep/MappedHumanGenesAverage2D' ids{idnr} 'id300'],'mappedX','color_RGB');
    datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' ids{idnr} '/'];
    load([datadir 'sample']);
    
    minx = min(mappedX);
    maxx = max(mappedX);
    maxval = max([abs(minx) abs(maxx)]);
    mappedab = mappedX*127/maxval;
    tC = makecform('lab2srgb');
    mappedRGB = applycform([50*ones(size(mappedab,1),1) mappedab],tC);
    
    % Plot figure
    nsamples = size(mappedX,1);
    figure(1)
    clf
    hold on
    for snr = 1:nsamples
        %     if ~isempty(colors.voxel.color_HEX{snr})
        %         colorRGB = rgbconv(colors.voxel.color_HEX{snr});
        %     else
        %         colorRGB = [0 0 0];
        %     end
        plot(mappedX(snr,1),mappedX(snr,2),'.','Color',mappedRGB(snr,:),'Markersize',markersize);
        %     text(mappedX(snr,1),mappedX(snr,2),sample.structure_name{snr});
    end
    
    % Format the plot
    axis off equal tight
    set(gcf,'Color','white')
    dpi = 300;
    figsizeinch = [7 7]/2.54;
    set(gcf,'Units','pixel','Position',[50 50 dpi*figsizeinch])
    saveas(1,['Figure5' figsubnums2(idnr) '.png']);
end


