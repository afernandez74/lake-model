function climate_20thCent=create20thCenturyClimate_Spinup_24(n_spinupyears)

load('../Data/climateCompilation20thCent_24.mat')

load('../Data/climate_avgYr.mat')

%% fill RH
%use model to predict rh
predicted_rh_uncert=normrnd(predicted_rh,(predicted_rh_ci(:,2)-predicted_rh_ci(:,1))/2);
predicted_rh_uncert(predicted_rh>100)=100;
rh_20thCent=[predicted_rh_uncert;agrimetcastor_rh(:,4)];
climate_20thCent=[final_pt,rh_20thCent];

%% fill Rs data
%create 20th century insolation (Rs) by looping measured data

Rs20thCentloopedwholes=readmatrix('../Data/Rs_20thCent_looped_w_holes_24.csv');
Rs=fillmissing(Rs20thCentloopedwholes(:,4),'linear');
climate_20thCent(:,7)=Rs;

%% fill Tw
%use model to predict Tw
predicted_tw_uncert=normrnd(predicted_tw,(predicted_tw_ci(:,2)-predicted_tw_ci(:,1))/2);
predicted_tw_uncert(predicted_tw<1)=1;

Tw_20thCent=[predicted_tw_uncert;castor_tw(:,4)];
climate_20thCent(:,8)=Tw_20thCent;

%% create spin-up time avg climate data 
if n_spinupyears~=0
    
    spinupyears=(climate_20thCent(1,3)-n_spinupyears : climate_20thCent(1,3)-1)';
    dates_spinup=(datetime(spinupyears(1),1,1):datetime(climate_20thCent(1,3),climate_20thCent(1,2),climate_20thCent(1,1))-1)';

    climate_spinup=NaN(size(dates_spinup,1),8);
    climate_spinup(:,1)=day(dates_spinup);
    climate_spinup(:,2)=month(dates_spinup);
    climate_spinup(:,3)=year(dates_spinup);
    j=1;

    for i = dates_spinup(1):dates_spinup(end)

        climate_spinup(j,4)=climate_avg(climate_avg(:,1)==day(i) & climate_avg(:,2)==month(i),4);
        climate_spinup(j,5)=climate_avg(climate_avg(:,1)==day(i) & climate_avg(:,2)==month(i),5);
        climate_spinup(j,6)=climate_avg(climate_avg(:,1)==day(i) & climate_avg(:,2)==month(i),6);
        climate_spinup(j,7)=climate_avg(climate_avg(:,1)==day(i) & climate_avg(:,2)==month(i),7);
        climate_spinup(j,8)=climate_avg(climate_avg(:,1)==day(i) & climate_avg(:,2)==month(i),8);
        j=j+1;
    end

    climate_20thCent=[climate_spinup;climate_20thCent];
end
end

