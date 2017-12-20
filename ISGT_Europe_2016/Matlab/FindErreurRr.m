function [ erreur_ru_rd ] = FindErreurRr(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ru_vec = linspace(0,2.5,26);
rd_vec = linspace(0,1,11);
erreur_ru_rd = NaN(length(ru_vec),length(rd_vec));


for i=1:1:length(ru_vec)
    ru = ru_vec(i);
    
    for j=1:1:length(rd_vec)
        rd = rd_vec(j);
        
        [Rr_matrix, Rr_max, Tr_opt, X_opt, Tr_theo, Rr_theo, X_theo] = RrResponseToTs( ru, rd );
        for ts_id=1:1:length(Rr_max)
            if Rr_max(ts_id)>Rr_theo(ts_id)
            erreur_ru_rd(i,j)=1

            end
        end

    end

end


end

