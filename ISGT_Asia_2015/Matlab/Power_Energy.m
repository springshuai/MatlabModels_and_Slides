function [ P_recharge, Bat_recharge, P_regulation, Bat_regulation ] = Power_Energy(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
P_d = 20;
P_n = 0.8*P_d;
rho_d = 0.45;
rho_u = rho_d;
rho_n = 1-rho_d-rho_u;
C_B = 50;
Delta = 0.1;

Num = 60;

Bat_recharge = zeros(1,Num);
Bat_regulation = zeros(1,Num);

for i = 1:1:Num
    
    if Bat_recharge(i)<C_B
        P_recharge(i) = P_d;
        Bat_recharge(i+1) = Bat_recharge(i)+P_recharge(i)*Delta;
        Bat_recharge(i+1) = min(Bat_recharge(i+1), C_B);
    else
        P_recharge(i) = 0;
        Bat_recharge(i+1) = Bat_recharge(i);
    end
    
    
    if Bat_regulation(i)<C_B
        signal = rand();
        if signal <= 0.45
            P_regulation(i) = 0;
            Bat_regulation(i+1) = Bat_regulation(i) + P_regulation(i)*Delta;
        elseif signal >= 1-0.45
            P_regulation(i) = P_d;
            Bat_regulation(i+1) = Bat_regulation(i) + P_regulation(i)*Delta;
            Bat_regulation(i+1) = min(Bat_regulation(i+1), C_B);
        else
            P_regulation(i) = P_n;
            Bat_regulation(i+1) = Bat_regulation(i) + P_regulation(i)*Delta; 
            Bat_regulation(i+1) = min(Bat_regulation(i+1), C_B);
        end
    else
            P_regulation(i) = 0;
            Bat_regulation(i+1) = Bat_regulation(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% -----write the data into TEX file Power_Energy.tex ----%%%%%%%
Parameters = ['$C_B = $' num2str(C_B) 'kWh, $P_d = $' num2str(P_d) 'kW, $P_n = $' num2str(P_n)...
    'kW, $\Delta = $' num2str(Delta) 'hour, $\rho_u=\rho_d=$' num2str(rho_u)];
x_label = 'Number of regulation sessions';
y_label = 'Charging power and energe obtained';
name_file = 'Power_Energy.tex';
fid = fopen(['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' name_file], 'w');

% fprintf(fid, '%s\n', Parameters);
fprintf(fid, '%s\n', ['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={1,1},anchor=south west},width=\figwidth,height=\figheight,cycle list name=\mylist,xlabel=' ...
    x_label ',ylabel=' y_label ']']);

fprintf(fid,'%s\n','\addplot +[const plot]coordinates{');
for i = 1:1:Num
    fprintf(fid, '%s', '(',num2str(i),',',num2str(P_recharge(i)),')');
end
fprintf (fid, '%s\n', '};');

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:Num
    fprintf(fid, '%s', '(',num2str(i),',',num2str(Bat_recharge(i)),')');
end
fprintf (fid, '%s\n', '};');

fprintf(fid,'%s\n','\addplot +[const plot]coordinates{');
for i = 1:1:Num
    fprintf(fid, '%s', '(',num2str(i),',',num2str(P_regulation(i)),')');
end
fprintf (fid, '%s\n', '};');

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:Num
    fprintf(fid, '%s', '(',num2str(i),',',num2str(Bat_regulation(i)),')');
end
fprintf (fid, '%s\n', '};');

fprintf (fid, '%s\n', '\end{axis}\end{tikzpicture}}');
fprintf (fid, '%s\n','\caption{Power and energy an EV obained with and without regulation ', Parameters, '}');
fclose (fid);
    
end

