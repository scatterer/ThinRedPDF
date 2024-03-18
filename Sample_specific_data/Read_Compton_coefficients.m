function [Compton_coeff] = Read_Compton_coefficients(Z)
a     = zeros(5, length(Z));
b     = zeros(5, length(Z));
c     = zeros(1, length(Z));


for k = 1:length(Z)
    [a(:,k), b(:,k), c(1,k)] = find_compton(Z(k));
end

Compton_coeff = struct('a', transpose(a), 'b', transpose(b), 'c', transpose(c));

    function [a, b, c] = find_compton(Z)
        a = zeros(5,1);
        b = zeros(5,1);
        c = 0;
        
        fid   = fopen( 'Form_Factor_and_Elemental_data/Compton_ScatFactor.dat');
        file  = textscan(fid,'%s','delimiter','\n');
        lines = file{1};
        fclose(fid);
        
        
        for i = 1:length(lines)
            if size(lines{i}) < 3
            else
                dummy_line = split(lines{i});
                if strcmp(dummy_line{1}, '#S')
                    
                    if str2num(dummy_line{2}) == Z
                        dummy_line = split(lines{i + 3});
                        for j = 1:5
                            a(j) = str2double(dummy_line{j});
                            b(j) = str2double(dummy_line{j + 6});
                        end
                        c = str2double(dummy_line{6});
                        break
                    end
                end
            end
            
        end
    end

end