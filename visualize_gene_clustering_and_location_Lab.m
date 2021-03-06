%% Visualize clustered data
clear variables
close all

ids = {'9861','10021','12876','14380','15496','15697'};
nbrains = length(ids);
fnaddition = 'Selection';

bnr = 2;
id = ids{bnr};
datadir = [LUMCDATADIR 'tSNE_ABA/processed/'];
load([datadir 'MappedGenes' fnaddition id]);
rawdatadir = [LUMCDATADIR 'tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
load([rawdatadir 'sample']);

%% Range of values
minx = min(mappedX);
maxx = max(mappedX);
maxval = max([abs(minx) abs(maxx)]);

mappedab = mappedX*127/maxval;
nsamples = size(mappedab,1);

%% Plot mapping with Lab colors
tC = makecform('lab2srgb');
mappedrgb = applycform([75*ones(nsamples,1) mappedab],tC);

figure(1)
clf
hold on
for snr = 1:nsamples
    plot(mappedab(snr,1),mappedab(snr,2),'.','Color',mappedrgb(snr,:));
end
axis equal

figure(2)
clf
hold on
for snr = 1:nsamples
    plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',mappedrgb(snr,:));
end
view(3)
axis equal

figure(3)
clf
hold on
for snr = 1:nsamples
    c = sample.color_RGB(snr,:);
    plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
end
axis equal

%%

% for bnr = 1:nbrains
% %% Load data
% id = ids{bnr};
% load(['MappedGenes' fnaddition id]);
% datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
% load([datadir 'sample']);
% 
% %% Plot with left and right brain separately
% % Colorspacel
% mincoords = [-80 -120 -80]-5;
% maxcoords = [80 80 80]+5;
% crange = maxcoords-mincoords;
% % crange = 2*max(abs(coords));
% % mincoords = -max(abs(coords));
% nsamples = size(X,1);
% 
% %
% figure(1)
% clf
% hold on
% for snr = 1:nsamples
%     c = (coords(snr,:)-mincoords)./crange;
%     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
% end
% saveas(1,['allgenes' fnaddition id]);
% 
% %% Plot with left and right together
% % Colorspace
% symcoords = coords;
% symcoords(:,1) = abs(symcoords(:,1));
% 
% %
% figure(2)
% clf
% hold on
% for snr = 1:nsamples
%     c = (symcoords(snr,:)-mincoords)./crange;
%     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
% end
% saveas(2,['allgenessymmetric' fnaddition id])
% 
% %% Show locations in brain
% figure(3)
% clf
% hold on
% for snr = 1:nsamples
%     c = (coords(snr,:)-mincoords)./crange;
%     plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',c);
% end
% view(-96,8)
% axis equal
% xlabel('x (MNI152)')
% ylabel('y (MNI152)')
% zlabel('z (MNI152)')
% saveas(3,['brainsamplesside' fnaddition id]);
% view(-180,0);
% saveas(3,['brainsamplesfront' fnaddition id]);
% 
% figure(4)
% clf
% hold on
% for snr = 1:nsamples
%     c = (symcoords(snr,:)-mincoords)./crange;
%     plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',c);
% end
% view(-96,8)
% axis equal
% xlabel('x (MNI152)')
% ylabel('y (MNI152)')
% zlabel('z (MNI152)')
% saveas(4,['brainsamplessidesymmetric' fnaddition id]);
% view(-180,0);
% saveas(4,['brainsamplesfrontsymmetric' fnaddition id]);
% 
% %% Region dependent coloring
% figure(5)
% clf
% hold on
% for snr = 1:nsamples
%     c = sample.color_RGB(snr,:);
%     plot(mappedX(snr,1),mappedX(snr,2),'.','Color',c);
% end
% saveas(5,['allgenesregions' fnaddition id]);
% 
% figure(6)
% clf
% hold on
% for snr = 1:nsamples
%     c = sample.color_RGB(snr,:);
%     plot3(coords(snr,1),coords(snr,2),coords(snr,3),'.','Color',c);
% end
% view(-96,8)
% axis equal
% xlabel('x (MNI152)')
% ylabel('y (MNI152)')
% zlabel('z (MNI152)')
% saveas(6,['brainsamplessideregions' fnaddition id]);
% view(-180,0);
% saveas(6,['brainsamplesfrontregions' fnaddition id]);
% 
% end

