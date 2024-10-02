function F_E = F_E (Ta,SA,E)

%Evaporation from lake surface m3/month

    if Ta<=0
        
        F_E=0;
    
    else
        
        F_E = SA * E;
    
    end
        
end

