function F_SF = F_SF (Ta,P,CA)
% Flux of snowfall over catchment

    if Ta > 0 
        
        F_SF = 0;
        
    else
        
        F_SF = (P/1000) * CA ;
    
    end
    

end

