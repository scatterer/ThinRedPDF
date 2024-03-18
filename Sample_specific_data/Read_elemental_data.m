function [Z, M] = Read_elemental_data(E)

Z = zeros(1, length(E));
M = zeros(1, length(E));

fid   = fopen( 'Form_Factor_and_Elemental_data/PERIODIC_TABLE.dat');
file  = textscan(fid,'%s','delimiter','\n');
lines = file{1};
fclose(fid);
lines = lines(29:end - 22);

for q = 1:length(E)
    for i = 1:length(lines)
        dummy_line = split(lines{i});
        
        if strcmp(dummy_line(3), E(q))
            Z(q) = str2double(dummy_line(1));
            M(q) = str2double(dummy_line(4));
        end
    end
end
end
