function res=for_wenjing(t_min,t_max,nb_points,freq_vec)
global string_command;
string_command=['for_wenjing(' num2str(t_min) ',' num2str(t_max) ',' int2str(nb_points) ',' num2str(freq_vec) ')'];


t_vec=linspace(t_min,t_max,nb_points);
for i=1:size(freq_vec,2)
	freq=freq_vec(i);
	for j=1:size(t_vec,2)
		t=t_vec(j);
		my_matrix(i,j)=sin(freq*t);
	end
end

ResultsToLaTex(my_matrix,t_vec,'legend entries={small period,larger period,even larger period,largest period}','time $t$','the result','my_filename.tex')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=ResultsToLaTex(Matrix,x_vector,axis_options,x_label,y_label,filename)
global string_command
fid=fopen(['results/' filename],'w');
c=clock;
infos=['%Obtained from ' string_command ', run on ' int2str(c(3)) '/' int2str(c(2)) '/' int2str(c(1)) ' at ' int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6))];
fprintf(fid,'%s\n',infos);
fprintf(fid,'%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={(1,1.03)},anchor=south east},width=\figwidth,height=\figheight,cycle list name=\mylist,every axis legend/.append style={nodes={right}},xlabel=' x_label ',ylabel=' y_label ',' axis_options ']']);

vectors=cell(1,size(Matrix,1));

for j=1:size(Matrix,1)
	for i=1:size(Matrix,2)
		vectors{j}=[vectors{j} '(' num2str(x_vector(i)) ',' num2str(Matrix(j,i)) ')'];
	end
	fprintf(fid,'%s\n',['\addplot coordinates{']);
	fprintf(fid,'%s\n',[vectors{j} '};']);
end

fprintf(fid,'%s\n','\end{axis}\end{tikzpicture}}');
fclose(fid);


