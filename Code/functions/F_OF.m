function F_OF = F_OF ( RES_SL, RES_DL,Max_V )

%Outflow flux m3/month

    if RES_SL+RES_DL > Max_V
        
        F_OF = (RES_SL+RES_DL-Max_V);
        
    else 
        
        F_OF=0;
    
    end
    
end
