function  F_DLM_daily  = F_DLM_daily (SVC,RES_SL,month)

% Deep lake mixing flux m3

    if (SVC-RES_SL) > 0 && month <= 11
        
        
        F_DLM_daily = SVC - RES_SL ;
    
    else
        
        F_DLM_daily=0;
    
    end

end