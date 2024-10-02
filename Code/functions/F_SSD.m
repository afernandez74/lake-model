function F_SSD = F_SSD (RES_SS,SSC)

% Surface soil drainage m3/month

    if RES_SS - SSC  > 0

        F_SSD = RES_SS - SSC;
    
    else
        
        F_SSD = 0;
    
    end
    
end
