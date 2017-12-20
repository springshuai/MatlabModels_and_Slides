function [ x_mono ] = RuRdPlane_mono(  )
%RURDPLANE Summary of this function goes here
%   Detailed explanation goes here

ru_mono_vec = linspace(0.0,2.5,251);
rd_mono_vec = linspace(0.0,1,101);
x_mono = zeros(length(rd_mono_vec),length(ru_mono_vec));


for id_ru = 1:1:length(ru_mono_vec)
    ru = ru_mono_vec(id_ru);
    for id_rd = 1:1:length(rd_mono_vec)
        rd = rd_mono_vec(id_rd);
        [x_mono(id_rd,id_ru),] = MainFunction_Mono(ru, rd);
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% write the data into HUB_rd_ru_mono.tex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_name = 'RuRdPlane_mono';
file_name = 'HUB_rd_ru_mono.tex';
file_path = '/Users/shuaiwenjing/Dropbox/Wenjing/Competition_model/Graphics/' ;
fid = fopen([file_path file_name],'w');

fprintf(fid,'%s\n',['% from' func_name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%------------- the region where x=1
fprintf(fid,'%s\n','% the boarder where x=1');
fprintf(fid,'%s\n',' coordinates{');

for id_rd = 1:1:length(rd_mono_vec)
    rd = rd_mono_vec(id_rd);
    for id_ru = 2:1:length(ru_mono_vec)
        ru = ru_mono_vec(id_ru);
        if (x_mono(id_rd,id_ru)==1)&&(x_mono(id_rd,id_ru-1)~=1)
            
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

for id_rd = 1:1:length(rd_mono_vec)
    rd = rd_mono_vec(id_rd);
    for id_ru = 1:1:length(ru_mono_vec)-1
        ru = ru_mono_vec(id_ru);
        if (x_mono(id_rd,id_ru)==0)&&(x_mono(id_rd,id_ru+1)~=0)
            
            fprintf(fid,'%s','(',num2str(ru),',',num2str(rd),')');
        end    
    end
end
fprintf(fid,'%s\n', '};');


fclose(fid);
end

