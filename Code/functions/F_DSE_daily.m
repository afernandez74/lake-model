function F_DSE_daily = F_DSE_daily(RES_DS,PET_demand,F_SSE,dt,TSC)

% Deep soil evapotranspiration flux m3/month
      if PET_demand * dt > F_SSE
          
          if RES_DS > PET_demand * dt - F_SSE
              
              F_DSE_daily = (PET_demand * dt - F_SSE)*RES_DS/TSC;
              
          else
              
              F_DSE_daily = RES_DS*dt;
              
          end
      else
          
          F_DSE_daily=0;
      end
end

