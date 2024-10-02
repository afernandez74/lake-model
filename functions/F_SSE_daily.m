function F_SSE_daily = F_SSE_daily (RES_SS,PET_demand,dt)
    
% Surface soil evapotranspiration flux m3/month
            
    if RES_SS > PET_demand * dt

        F_SSE_daily = PET_demand * dt ;

    else

        F_SSE_daily = RES_SS * dt ;

    end

end

