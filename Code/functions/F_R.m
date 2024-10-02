function F_R = F_R (Ta,P,CA)

% Rainfall over catchment m3

    if Ta <= 0 
        
        F_R = 0;
        
    else
        
        F_R = (P/1000) * CA ;
    
    end
    
end

