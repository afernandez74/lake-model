function  F_SLM_daily  = F_SLM_daily (RES_SL,month)

%Surface lake mixing flux m3/month

    if month < 10
        
        F_SLM_daily = 0;
    
    else
        
        F_SLM_daily = RES_SL ;
    
    end
    
end

