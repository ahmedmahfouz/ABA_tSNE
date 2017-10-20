%% Align ABA, based on samples to ABA
clear variables
close all

% Add Nifti reader
addpath('../NIFTI/');
addpath('../tracking/MetaIO/');

% Load ABA data
rawdatadir = [LUMCDATADIR 'tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor9861/'];
sampledata = load([rawdatadir 'sample']);

% Load ABA brain
brainfile = [rawdatadir 'T1.nii'];
brain = load_nii(brainfile);

% Load brain atlas
atlasfile = [LUMCDATADIR '/tSNE_ABA/MNI_atlas/icbm_avg_152_t1_tal_nlin_symmetric_VI.nii'];
atlas = load_nii(atlasfile);
atlasmask = load_nii([LUMCDATADIR '/tSNE_ABA/MNI_atlas/icbm_avg_152_t1_tal_nlin_symmetric_VI_mask.nii']);
atlas.img(atlasmask.img==0) = min(atlas.img(:));

% %% Convert ABA samples to brain atlas
% coords = [sampledata.sample.mni_x,...
%     sampledata.sample.mni_y,...
%     sampledata.sample.mni_z];
% A = [atlas.hdr.hist.srow_x;
%     atlas.hdr.hist.srow_y;
%     atlas.hdr.hist.srow_z];
% Ahom = [A; 0 0 0 1];
% imgcoords = Ahom\[coords ones(size(coords,1),1)]';
% imgcoords = imgcoords(1:3,:)'+1;
% 
% abaimg = false(size(atlas.img));
% imgind = sub2ind(size(abaimg),round(imgcoords(:,1)),round(imgcoords(:,2)),round(imgcoords(:,3)));
% abaimg(imgind) = true;
% abaimg = bdilation(abaimg,10);

%% Register both brains using elastix
scaledatlas = (atlas.img-min(atlas.img(:)))/(max(atlas.img(:))-min(atlas.img(:)));
metaImageWrite(scaledatlas,'atlas.mhd');
metaImageWrite(brain.img,'brain.mhd');

% Register atlas to brain
elastixcmd = '"C:\Program Files (x86)\elastix_v4.7\elastix"';
transformixcmd = '"C:\Program Files (x86)\elastix_v4.7\transformix"';
regcmd = [elastixcmd sprintf(' -f %s -m %s -out . -p Par0000affine.txt -p Par0000bspline.txt','brain.mhd','atlas.mhd')];
system(regcmd);

%% Transport coordinates to MNI brain
trcmd = [transformixcmd sprintf(' -in %s -out . -tp %s','atlas.mhd','TransformParameters.1.txt')];
system(trcmd);
regimg = metaImageRead('result.mhd');

%% Get new coordinates in MNI space (voxels in ABA brain are 1mm isotropic)
imgcoords = double([size(brain.img,3)-sampledata.sample.mri_voxel_z,...
    size(brain.img,2)-sampledata.sample.mri_voxel_x,...
    size(brain.img,1)-sampledata.sample.mri_voxel_y]);
imgcoords = imgcoords+repmat([12 -7 -7],[size(imgcoords,1) 1]);
% imgcoords = [sampledata.sample.mri_voxel_x,...
%     sampledata.sample.mri_voxel_y,...
%     sampledata.sample.mri_voxel_z];
% Ab = [brain.hdr.hist.srow_x;
%     brain.hdr.hist.srow_y;
%     brain.hdr.hist.srow_z];
% Abhom = [Ab; 0 0 0 1];
% c = Abhom*[double(imgcoords'); ones(1,size(imgcoords,1))];

write_elastix_points(imgcoords,'expimgcoords.txt');

trptscmd = [transformixcmd sprintf(' -def %s -out . -tp %s','expimgcoords.txt','TransformParameters.1.txt')];
system(trptscmd);

[outputpoints, inputpoints] = read_elastix_points('outputpoints.txt');

%% Plot transformed points
figure(1)
clf
hold on
plot3(outputpoints(:,1),outputpoints(:,2),outputpoints(:,3),'.');
plot3(inputpoints(:,1),inputpoints(:,2),inputpoints(:,3),'r.');
view(3)
axis equal

%% Project points back into the image
icind = sub2ind(size(brain.img),imgcoords(:,2),...
    imgcoords(:,1),...
    imgcoords(:,3));

icimg = false(size(brain.img));
icimg(icind) = true;

%% Project points into the atlas
opind = sub2ind(size(atlas.img),round(outputpoints(:,2)),...
    round(outputpoints(:,1)),...
    round(outputpoints(:,3)));

opimg = false(size(atlas.img));
opimg(opind) = true;
