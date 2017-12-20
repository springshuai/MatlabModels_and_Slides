function [R_compare, R_matrix ] = Print_R_matrix(  )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
gamma = 0.05;
r_u = 2;
r_d = 0.6;
rho_u = 0.48;
rho_d = 0.48;

[ R_compare, U_bar, w_P_vec, Ar_vec, Ac_vec, R_matrix, P_A] = R_Tc_Tr_Correction( gamma, r_u, r_d, rho_u, rho_d );

file_ccv = 'R_3D_raw_ccv.tex';
file_non = 'R_3D_raw_non.tex';
fid_ccv = fopen(['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' file_ccv], 'w');
fid_non = fopen(['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' file_non], 'w');

fprintf (fid_ccv, '%s\n', '% gamma=0.05, r_u=2.0, r_d=0.6, rho_u=0.48, rho_d=0.48 )');

Num = size(R_matrix);
for i = 1:1:Num
    for j=1:1:Num
        if R_matrix(i,j,1)/20>R_matrix(i,j,2)/P_A
            fprintf(fid_ccv, '%s','(',num2str(R_matrix(i,j,1)),',', num2str(R_matrix(i,j,2)),',', num2str(R_matrix(i,j,3)),')');
            fprintf(fid_non, '%s','(',num2str(R_matrix(i,j,1)),',', num2str(R_matrix(i,j,2)),',', num2str(NaN), ')');
        else
            fprintf(fid_ccv, '%s','(',num2str(R_matrix(i,j,1)),',', num2str(R_matrix(i,j,2)),',', num2str(NaN),')');
            fprintf(fid_non, '%s','(',num2str(R_matrix(i,j,1)),',', num2str(R_matrix(i,j,2)),',', num2str(R_matrix(i,j,3)), ')');
         end
    end
end
fclose(fid_ccv);
fclose(fid_non);
        

end

