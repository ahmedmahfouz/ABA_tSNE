%% Write a HTML file with the region colors
clear variables
close all

id = '9861';
datadir = ['/home/mvandegiessen/data/tSNE_ABA/rawData_25Feb2014/normalized_microarray_donor' id '/'];
load([datadir 'sample']);

% Get list of unique structures
[C,IC] = unique(sample.structure_id);
structure_id = C;
structure_HEX = sample.color_HEX(IC);
structure_name = sample.structure_name(IC);


% Output header
fid = fopen('region_colors.html','w');
fprintf(fid,'<!DOCTYPE html>\n<html>\n<body>\n');
fprintf(fid,'<table>');

for snr = 1:length(C)
    fprintf(fid,'<tr>\n');
    fprintf(fid,'<td>%d</td><td>%s</td><td bgcolor="#%s"><font color="%s">XXXX</font></td>\n',...
        structure_id(snr),structure_name{snr},structure_HEX{snr},structure_HEX{snr});
    fprintf(fid,'</tr>\n');
end

% Output footer
fprintf(fid,'</table>');
fprintf(fid,'</body>\n</html>\n');
fclose(fid);