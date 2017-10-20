%% Map gene tsne colors to brodmann areas
clear variables
close all

% Add Nifti reader
addpath('..\NIFTI');

% Load atlas areas
% atlas = load_nii([LUMCDATADIR '/mricron/templates/brodmann.nii']);
% atlas = load_nii([LUMCDATADIR '/mricron/templates/aal.nii']);
atlas = load_nii([LUMCDATADIR '/mricron/templates/ch2bet.nii']);
atlas.img = (atlas.img>40)*1.0;
% atlas = load_nii([LUMCDATADIR '/mricron/templates/JHU-WhiteMatter-labels-1mm.nii']);
nlabels = max(atlas.img(:));

A = [atlas.hdr.hist.srow_x;
    atlas.hdr.hist.srow_y;
    atlas.hdr.hist.srow_z];
Ahom = [A; 0 0 0 1];

coord = [0 0 0];
imgcoord = Ahom\[coord 1]'+[1 1 1 0]';

% Load gene expression data
id = '9861';
% id = '10021';
procdatadir = [LUMCDATADIR 'tSNE_ABA/normalized/results_human_average_300/'];
mappeddata = load([procdatadir 'MappedHumanGenesAverage2D' id 'id300']);
rawdatadir = [LUMCDATADIR 'tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
sampledata = load([rawdatadir 'sample']);

%% Project sample coordinates to atlas areas
expcoords = [sampledata.sample.mni_x,...
    sampledata.sample.mni_y,...
    sampledata.sample.mni_z];
expvals = mappeddata.mappedX;
% expvals = (expvals-repmat(min(expvals),[size(expvals,1) 1]))./repmat(max(expvals)-min(expvals),[size(expvals,1) 1]);
% expimgcoord = Ahom\[expcoords ones(size(expcoords,1),1)]'+repmat([1 1 1 0]',[1 size(expcoords,1)]);

manualadjust = [0 20 -20];

isz = size(atlas.img);
expimgcoord = expcoords+repmat((isz-1)/2+manualadjust,[size(expcoords,1) 1])+1;
expimgcoord(expimgcoord<1) = 1;
expimgind = sub2ind(isz,round(expimgcoord(:,1)),round(expimgcoord(:,2)),round(expimgcoord(:,3)));

ovimg = zeros(isz);
ovimg(expimgind) = 1;

G = griddedInterpolant(double(atlas.img),'nearest','none');
explabs = G(expimgcoord);

%% Interpolant for full brain
ndims = size(expvals,2);
labind = find(atlas.img>0);
[labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
labicoords = [labix,labiy,labiz];
% labcoords = (A*[labicoords'; ones(1,size(labicoords,1))])';

fimg = cell(1,ndims);
for dnr = 1:ndims
    fimg{dnr} = zeros(size(atlas.img));
end
labvals = zeros(size(labicoords,1),ndims);
for dnr = 1:ndims
    F = scatteredInterpolant(expimgcoord,expvals(:,dnr),'nearest','nearest');
    labvals(:,dnr) = F(labicoords);
    fimg{dnr}(labind) = labvals(:,dnr);
end

%% Map to correct colors
minx = min(expvals);
maxx = max(expvals);
maxval = max([abs(minx) abs(maxx)]);
mappedab = [fimg{1}(labind) fimg{2}(labind)]*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([50*ones(size(mappedab,1),1) mappedab],tC);

colorimgr = zeros(size(atlas.img));
colorimgg = zeros(size(atlas.img));
colorimgb = zeros(size(atlas.img));
colorimgr(labind) = mappedRGB(:,1);
colorimgg(labind) = mappedRGB(:,2);
colorimgb(labind) = mappedRGB(:,3);
colorimg = joinchannels('rgb',colorimgr,colorimgg,colorimgb);

%% Interpolant for region based
maxlab = max(atlas.img(:));

limg = cell(1,ndims);
for dnr = 1:ndims
    limg{dnr} = zeros(size(atlas.img));
end

for lnr = 1:maxlab
    lnr
    % Get coordinates within the lab
    labind = find(atlas.img==lnr);
    [labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
    labicoords = [labix,labiy,labiz];
    labvals = zeros(size(labicoords,1),ndims);
    
    % Get coordinates that belong to this label
    mapsel = explabs==lnr;
    
    if sum(mapsel)>0
        pts = expimgcoord(mapsel,:);
        if size(pts,1) > 3
            for dnr = 1:ndims
                ev = expvals(mapsel,dnr);
                F = scatteredInterpolant(pts,ev,'linear','nearest');
                labvals(:,dnr) = F(labicoords);
                limg{dnr}(labind) = labvals(:,dnr);
            end
        else
            npt = dsearchn(pts,labicoords);
            ev = expvals(mapsel,:);
%             D = delaunayTriangulation(pts);
%             npt = nearestNeighbor(D,labicoords);
            for dnr = 1:ndims
                limg{dnr}(labind) = ev(npt,dnr);
            end
        end
    end
end

%% Map to correct colors
labind = find(atlas.img>0);
minx = min(expvals);
maxx = max(expvals);
maxval = max([abs(minx) abs(maxx)]);
limg{1}(isnan(limg{1})) = 0;
limg{2}(isnan(limg{2})) = 0;
mappedab = [limg{1}(labind) limg{2}(labind)]*127/maxval;
tC = makecform('lab2srgb');
mappedRGB = applycform([75*ones(size(mappedab,1),1) mappedab],tC);

colorimgr = zeros(size(atlas.img));
colorimgg = zeros(size(atlas.img));
colorimgb = zeros(size(atlas.img));
colorimgr(labind) = mappedRGB(:,1);
colorimgg(labind) = mappedRGB(:,2);
colorimgb(labind) = mappedRGB(:,3);
lcolorimg = joinchannels('rgb',colorimgr,colorimgg,colorimgb);

%% Interpolant for full brain with original colors
ndims = size(mappeddata.color_RGB,2);
labind = find(atlas.img>0);
[labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
labicoords = [labix,labiy,labiz];
% labcoords = (A*[labicoords'; ones(1,size(labicoords,1))])';

rimg = cell(1,ndims);
for dnr = 1:ndims
    rimg{dnr} = zeros(size(atlas.img));
end
labvals = zeros(size(labicoords,1),ndims);
for dnr = 1:ndims
    F = scatteredInterpolant(expimgcoord,mappeddata.color_RGB(:,dnr),'nearest','nearest');
    labvals(:,dnr) = F(labicoords);
    rimg{dnr}(labind) = labvals(:,dnr);
end

rcolorimg = joinchannels('rgb',rimg{1},rimg{2},rimg{3});

break

%% Select coordinates from a label
nvals = size(expvals,2);
pimg = cell(1,nvals);
for vnr = 1:nvals
    pimg{vnr} = zeros(size(atlas.img));
end

for lnr = 1:nlabels
    % Select coordinates from label
    sellab = lnr;
    labind = find(atlas.img==sellab);
    [labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
    labicoords = [labix-1 labiy-1 labiz-1];
    
    % Transfer to MNI space
    labcoords = (A*[labicoords'; ones(1,size(labicoords,1))])';
    
    % Select genecoordinates with the same label
    selexpcoords = expcoords(explabs==sellab,:);
    selexpvals = expvals(explabs==sellab,:);
    
    % Interpolate between datapoints
    if size(selexpcoords,1) > 0
%         labvals = zeros(size(labcoords,1),nvals);
%         for dnr = 1:nvals;
%             F = scatteredInterpolant(selexpcoords,selexpvals(:,dnr),'nearest');
%             labvals(:,dnr) = F(labcoords);
%             pimg{dnr}(labind) = mean(labvals(:,dnr));
%         end
        for dnr = 1:nvals
            pimg{dnr}(labind) = median(selexpvals(:,dnr));
        end
    end
end

cpimg = joinchannels('rgb',pimg{1},pimg{2},pimg{3});
labmask = atlas.img>0;
lpimg = joinchannels('lab',labmask*1.0,pimg{2},pimg{1});

%% Check coordinate spaces
% figure(10)
% clf
% labind = find(atlas.img>0);
% [labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
% labicoords = [labix-1 labiy-1 labiz-1];
% labcoords = (A*[labicoords'; ones(1,size(labicoords,1))])';
% plot3(labcoords(:,1),labcoords(:,2),labcoords(:,3),'.');
% view(3)
% axis equal
% hold on
% plot3(expcoords(:,1),expcoords(:,2),expcoords(:,3),'r.');

%% Get interpolated results for the whole brain
labind = find(atlas.img>0);
[labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
labicoords = [labix-1 labiy-1 labiz-1];
labcoords = (A*[labicoords'; ones(1,size(labicoords,1))])';

fimg = cell(1,nvals);
for vnr = 1:nvals
    fimg{vnr} = zeros(size(atlas.img));
end
labvals = zeros(size(labcoords,1),nvals);
for dnr = 1:nvals
    F = scatteredInterpolant(expcoords,expvals(:,dnr),'nearest','nearest');
    labvals(:,dnr) = F(labcoords);
    fimg{dnr}(labind) = labvals(:,dnr);
end

%% Show image
fcimg = joinchannels('rgb',fimg{1},fimg{2},fimg{3});

%% Project on brain surface
% Select half brain
labmaskhalf = labmask;
labmaskhalf(ceil(size(labmaskhalf,1)/2):end,:,:) = false;
[faces,ivertices] = isosurface(gaussf(labmaskhalf,1));
ivertices = [ivertices(:,2) ivertices(:,1) ivertices(:,3)];
vertices = (A*[ivertices';ones(1,size(ivertices,1))])';

vvals = zeros(size(vertices,1),nvals);
for dnr = 1:nvals
    F = scatteredInterpolant(expcoords,expvals(:,dnr),'nearest','nearest');
    vvals(:,dnr) = F(vertices);
end

%% Show vertices
figure(10)
clf
patch('Faces',faces,'Vertices',vertices,'FaceColor','interp',...
    'EdgeColor','none','FaceVertexCData',[vvals(:,1) vvals(:,2) vvals(:,3)]);
view(3)
axis equal
camlight
material dull

% %% Image with sample colors from ABA
% labind = find(atlas.img>0);
% [labix,labiy,labiz] = ind2sub(size(atlas.img),labind);
% labicoords = [labix-1 labiy-1 labiz-1];
% labcoords = (A*[labicoords'; ones(1,size(labicoords,1))])';
% 
% nlvals = size(sampledata.sample.color_RGB,2);
% 
% abalimg = cell(1,nlvals);
% for vnr = 1:nlvals
%     abalimg{vnr} = zeros(size(atlas.img));
% end
% labvals = zeros(size(labcoords,1),nlvals);
% for dnr = 1:nlvals
%     F = scatteredInterpolant(expcoords,sampledata.sample.color_RGB(:,dnr),'nearest','nearest');
%     labvals(:,dnr) = F(labcoords);
%     abalimg{dnr}(labind) = labvals(:,dnr);
% end
% 
% %% Show image
% abalimg = joinchannels('rgb',abalimg{1},abalimg{2},abalimg{3});
% 
% %% Project on brain surface
% [faces,ivertices] = isosurface(gaussf(labmask,1));
% ivertices = [ivertices(:,2) ivertices(:,1) ivertices(:,3)];
% vertices = (A*[ivertices';ones(1,size(ivertices,1))])';
% 
% vabalvals = zeros(size(vertices,1),nlvals);
% for dnr = 1:nlvals
%     F = scatteredInterpolant(expcoords,sampledata.sample.color_RGB(:,dnr),'nearest','nearest');
%     vabalvals(:,dnr) = F(vertices);
% end
% 
% %% Show vertices
% figure(11)
% clf
% patch('Faces',faces,'Vertices',vertices,'FaceColor','interp',...
%     'EdgeColor','none','FaceVertexCData',[vabalvals(:,1) vabalvals(:,2) vabalvals(:,3)]);
% view(3)
% axis equal
% camlight
% material dull
% lighting phong

