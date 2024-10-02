function F_SM_daily = F_SM_daily (Ta,CA,RES_SP,dt)

% Snow melt flux m3/month

        if Ta < -2

            F_SM_daily = 0;   

        elseif RES_SP > 0.7/1000 * (2+Ta) * CA

            F_SM_daily = 0.7/1000 * (2+Ta) * CA * dt;

        else 
            
            F_SM_daily = RES_SP * dt ;
            
        
        end
    
    
end

