%WRITE_ELASTIX_POINTS  Write elastix points to file

function result = write_elastix_points(pts,filename)

% Number of points
npts = size(pts,1);

% Open file
fid = fopen(filename,'w');
fprintf(fid,'point\n');
fprintf(fid,'%d\n',npts);
for pnr = 1:npts
    fprintf(fid,'%1.10f %1.10f %1.10f\n',pts(pnr,1),pts(pnr,2),pts(pnr,3));
end
fclose(fid);

result = 0;
