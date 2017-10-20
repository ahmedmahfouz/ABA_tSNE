function [outputpoints, inputpoints] = read_elastix_points(filename)

% Start with empty lists
outputpoints = [];
inputpoints = [];

fid = fopen(filename);
while ~feof(fid)
    l = fgetl(fid);
    s = textscan(l,'%s','Delimiter',';');
    ips = s{1}{3};
    ops = s{1}{5};
    outputpoints = [outputpoints;
        sscanf(ops,'OutputPoint = [ %f %f %f ]')'];
    inputpoints = [inputpoints;
        sscanf(ips,'InputPoint = [ %f %f %f ]')'];
end

fclose(fid);