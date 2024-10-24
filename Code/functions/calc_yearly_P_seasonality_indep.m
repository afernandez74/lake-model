%total precip for every year except for imcomplete years (first and last)
function yearly_P=calc_yearly_P_seasonality_indep(P_mod_summer,P_mod_winter,climate)
    
    year_o=climate(1,3);
    year=NaN(climate(end,3)-climate(1,3)-1,1);
    P_total=NaN(climate(end,3)-climate(1,3)-1,1);
    P_summer=NaN(climate(end,3)-climate(1,3)-1,1);
    P_winter=NaN(climate(end,3)-climate(1,3)-1,1);
    P_daily_change_summer=NaN(climate(end,3)-climate(1,3)-1,1);
    P_daily_change_winter=NaN(climate(end,3)-climate(1,3)-1,1);
    
    for i=1:(climate(end,3)-climate(1,3)-1)
        year_i=i+year_o;
        year(i)=year_i;
        P_total_year_i=sum(climate(climate(:,3)==(year_i),4));
        P_summer_i=sum(climate(climate(:,3)==(year_i) & climate(:,2)>3 & climate(:,2)<10,4));
        P_winter_i=P_total_year_i-P_summer_i;
        P_total(i)=P_total_year_i;
        P_summer(i)=P_summer_i;
        P_winter(i)=P_winter_i;
        
        P_daily_change_summer_i=P_mod_summer*(P_total_year_i/P_summer_i);
        P_daily_change_winter_i=P_mod_winter*(P_total_year_i/P_winter_i);
        P_daily_change_summer(i)=P_daily_change_summer_i;
        P_daily_change_winter(i)=P_daily_change_winter_i;
    end
    yearly_P=[year P_total P_summer P_winter P_daily_change_summer P_daily_change_winter];

end
