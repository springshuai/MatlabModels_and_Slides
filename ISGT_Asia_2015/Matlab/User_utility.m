function [  ] = User_utility( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Start = 0;
End = 0.6;
Step = 0.01;
Num = length(Start:Step:End);
c_b = 50;
p_d = 20;
p_a = 8;

% first group of parameters fall in the first case
t_c1 = 0.15;
t_r1 = 0.04;

% since t_r1/p_a <= t_c1/p_d
    for i = 1:1:Num
        theta(i) = Start + (i-1)*Step;
        if  theta(i) <= t_r1*c_b/p_a
            U1(i) = 0;
        elseif  theta(i) <= (t_c1-t_r1)*c_b/(p_d-p_a)
            U1(i) = theta(i)*p_a - t_r1*c_b;
        else 
            U1(i) = theta(i)*p_d - t_c1*c_b;
        end
    end
    
% second group of parameters fall in the second case
t_c2 = 0.15;
t_r2 = 0.10;
% since t_r2/p_a > t_c2/p_d
    for j = 1:1:Num
        theta(j) = Start + (j-1)*Step;
        if theta(j) > t_c2*c_b/p_d
            U2(j) = theta(j)*p_d - t_c2*c_b;
        else
            U2(j) = 0;
        end
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% -----write the data into TEX file User_utility.tex ----%%%%%%%

Parameters = ['$C_B = $' num2str(c_b) ', $P_d = $' num2str(p_d) ', $P_A = $'...
    num2str(p_a) ', $T_c = $' num2str(t_c1) ', $T_{r} = $' num2str(t_r1) ' or ' num2str(t_r2) ...
    ', $\theta$ goes from ' num2str(Start) ' to ' num2str(End)];
x_label = 'User preference parameter $\theta$';
y_label = 'User utility';
name_file = 'User_utility.tex';
fid = fopen(['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' name_file], 'w');

% fprintf(fid, '%s\n', Parameters);
fprintf(fid, '%s\n', ['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={1,1},anchor=south west},width=\figwidth,height=\figheight,cycle list name=\mylist,xlabel=' ...
    x_label ',ylabel=' y_label ']']);

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:Num
    fprintf(fid, '%s', '(',num2str(theta(i)),',',num2str(U1(i)),')');
end
fprintf (fid, '%s\n', '};');

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:Num
    fprintf(fid, '%s', '(',num2str(theta(i)),',',num2str(U2(i)),')');
end
fprintf (fid, '%s\n', '};');

fprintf (fid, '%s\n', '\end{axis}\end{tikzpicture}}');
fprintf (fid, '%s\n','\caption{User utility when ', Parameters, '}');
fclose (fid);

end

