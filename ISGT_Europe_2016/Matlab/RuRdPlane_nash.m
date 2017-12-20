function [ x_nash ] = RuRdPlane_nash(  )
%RURDPLANE_NASH Summary of this function goes here
%   Detailed explanation goes here

% rd_nash_vec = linspace(0.425,0.425,1);
ru_nash_vec = linspace(0,2.5,251);
rd_nash_vec = linspace(0,1,101);

x_nash = zeros(length(rd_nash_vec),length(ru_nash_vec));

for id_rd = 1:1:length(rd_nash_vec)
    rd = rd_nash_vec(id_rd);
    for id_ru = 1:1:length(ru_nash_vec)
%         ru = 1.1*rd+0.94+0.001*id_ru;
        ru = ru_nash_vec(id_ru);
        [x_nash(id_rd,id_ru)] = MainFunction_Nash(ru, rd);
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% write the data into HUB_rd_ru_nash.tex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_name = 'RuRdPlane_nash';
file_name = 'HUB_rd_ru_nash.tex';
file_path = '/Users/shuaiwenjing/Dropbox/Wenjing/Competition_model/Graphics/' ;
fid = fopen([file_path file_name],'w');

fprintf(fid,'%s\n',['% from' func_name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%------------- the region where x=1
fprintf(fid,'%s\n','% the boarder where x=1');
fprintf(fid,'%s\n',' coordinates{');

for id_rd = 1:1:length(rd_nash_vec)
    rd = rd_nash_vec(id_rd);
    for id_ru = 2:1:length(ru_nash_vec)
        ru = ru_nash_vec(id_ru);
        if (x_nash(id_rd,id_ru)==1)&&(x_nash(id_rd,id_ru-1)~=1)
            
            fprintf(fid,'%s','(',num2str(ru),',',num2str(rd),')');
        end    
    end
end
fprintf(fid,'%s\n', '};');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%------------- the region where x=0
fprintf(fid,'%s\n','% the boarder where x=0');
fprintf(fid,'%s\n',' coordinates{');

for id_rd = 1:1:length(rd_nash_vec)
    rd = rd_nash_vec(id_rd);
    for id_ru = 1:1:length(ru_nash_vec)-1
        ru = ru_nash_vec(id_ru);
        if (x_nash(id_rd,id_ru)==0)&&(x_nash(id_rd,id_ru+1)~=0)
            
            fprintf(fid,'%s','(',num2str(ru),',',num2str(rd),')');
        end    
    end
end
fprintf(fid,'%s\n', '};');


fclose(fid);
end

