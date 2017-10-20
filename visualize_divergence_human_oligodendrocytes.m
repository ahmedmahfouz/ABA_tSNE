% Compute overlap of samples per organ by computing the Jensen-Shannon
% divergence
clear variables
close all

%% Data
resultsdatadir = [LUMCDATADIR 'tSNE_ABA/normalized/results_human_average_cell_selection/'];
ids = {'9861','10021','12876','14380','15496','15697'};
id = ids{2};
datadir = [LUMCDATADIR 'tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor'];
% Load mapped data
load([resultsdatadir 'MappedHumanGenesAverage2Doligodendrocytes' id 'id92'],'mappedX');
% Load voxel labels
mainStructures = {'FL', 'PL', 'TL', 'OL', 'HiF', 'Str', 'GP', 'Amg', 'TH', ...
    'Hy', 'MES', 'Pons', 'MY', 'Cb', 'WM'};
voxelLabels = [];
load([datadir id '/sampleLables']);
voxelLabels = sampleLables(:);

nlabels = max(voxelLabels);

%% Determine histogram size
mappedX = mappedX(:,1:2);
minbin = min(mappedX);
maxbin = max(mappedX);
nbins = 20;

%% Compute divergence
D = zeros(nlabels,nlabels);

for lnr1 = 1:nlabels
    for lnr2 = 1:nlabels
        if lnr2 >= lnr1-1
            P = mappedX(voxelLabels==lnr1,:);
            Q = mappedX(voxelLabels==lnr2,:);
            d = jsdiv_pts(P,Q,nbins,minbin,maxbin);
            D(lnr1,lnr2) = d;
        else
            D(lnr1,lnr2) = NaN;
%             D(lnr1,lnr2) = 0;
        end
    end
end

binwidth = (maxbin-minbin)./nbins;
% save([resultsdatadir 'jsdiv_mouse_coronal_SVD'],'D','minbin','maxbin','nbins','mainStructures');

%% Show divergencies
Ddisp = [D zeros(size(D,1),1); nan(1,size(D,2)-1) 0 0];
fs = 15;
figure(1)
clf
hsurf = surface(Ddisp);
set(gcf,'Renderer','Painters')
set(gca,'FontSize',fs)
% set(gca,'XTick',1:nlabels,'XTickLabel',mainStructures);
% set(gca,'YTick',1:nlabels,'YTickLabel',mainStructures);
axis equal off
% set(gca,'YLim',[0.5 nlabels+0.5])
% rotateXLabels(gca,45)
cm = hot;
colormap(cm);
set(gcf,'Color','White')
rotate(hsurf,[0 0 1],-45)

textxloc = ((1:size(D,1))-size(D,1)/2)*sqrt(2)+size(D,1)/2-sqrt(2)/2+1;
textyloc = zeros(1,size(D,1))+size(D,1)/2+1.2+sqrt(2)/2;

for lnr = 1:size(D,1)
    text(textxloc(lnr),textyloc(lnr),mainStructures{lnr},'FontSize',fs,'Rotation',90)
end

set(gca,'CLim',[-0.01 max(D(:))])
h = colorbar('horiz');
set(h,'FontSize',fs)
set(h,'XTick',[0 10 20 30]);


set(gca,'XLim',[textxloc(1)-sqrt(2) textxloc(end)+sqrt(2)])
set(gca,'YLim',[textxloc(1)-sqrt(2) textyloc(end)+4])

set(gca,'CLim',[0 33]);

% orient landscape

set(gcf, 'InvertHardcopy', 'off','PaperPositionMode','auto');
saveas(1,'Figure9f.eps','epsc2');

% tightInset = get(gca, 'TightInset');
% position(1) = tightInset(1);
% position(2) = tightInset(2);
% position(3) = 1 - tightInset(1) - tightInset(3);
% position(4) = 1 - tightInset(2) - tightInset(4);
% set(gca, 'Position', position);

% %% Plot clusters
% lnr1 = 1;
% lnr2 = 5;
% 
% figure(2)
% clf
% hold on
% plot(mappedX(voxelLabels==lnr1,1),mappedX(voxelLabels==lnr1,2),'.');
% plot(mappedX(voxelLabels==lnr2,1),mappedX(voxelLabels==lnr2,2),'r.');
% lx = minbin(1):binwidth(1):maxbin(1);
% ly = minbin(2):binwidth(2):maxbin(2);
% plot([lx;lx],[ly(1)*ones(size(lx));ly(end)*ones(size(lx))],'k')
% plot([lx(1)*ones(size(ly));lx(end)*ones(size(ly))],[ly;ly],'k')
% 
% %%
% figure(3)
% clf
% hold on
% sl1 = find(voxelLabels==lnr1);
% sl2 = find(voxelLabels==lnr2);
% for snr = sl1
%     if ~isempty(voxel.color_HEX{snr})
%         c = rgbconv(voxel.color_HEX{snr});
%     else
%         c = [0 0 0];
%     end
%     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
% end
% for snr = sl2
%     if ~isempty(voxel.color_HEX{snr})
%         c = rgbconv(voxel.color_HEX{snr});
%     else
%         c = [0 0 0];
%     end
%     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
% end
% 
% plot([lx;lx],[ly(1)*ones(size(lx));ly(end)*ones(size(lx))],'k')
% plot([lx(1)*ones(size(ly));lx(end)*ones(size(ly))],[ly;ly],'k')
