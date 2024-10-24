function climate=create20thCenturyClimate_Spinup_20thCent_RandYrBlocks_24(n_spinupyears,n_reps)

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

climate_20thCent=climate_20thCent(climate_20thCent(:,3) > climate_20thCent(1,3) &...
    climate_20thCent(:,3) < max(climate_20thCent(:,3)),:);

%% re-order climate dataset in 10yr blocks
blocks=[1900 1909; 1910 1919; 1920 1929;1930 1939;1940 1949;1950 1959;...
    1960 1969;1970 1979;1980 1989;1990 1999;2000 2009;2010 2021];

blocks_rand = blocks(randperm(size(blocks,1)),:);
blocks_rand = repmat(blocks_rand,n_reps,1);
temp=NaN(1,8);

for i=1:length(blocks_rand)
    block_i_lo=blocks_rand(i,1);
    block_i_hi=blocks_rand(i,2);
    clim_block=climate_20thCent(climate_20thCent(:,3) >= block_i_lo ...
        & climate_20thCent(:,3) <= block_i_hi,:);
    temp=[temp;clim_block];
end

temp(1,:)=[];
clim_blocks=temp;
clim_blocks(1,3)=climate_20thCent(end,3)+1;
yr_change=find(clim_blocks(:,1)==1 & clim_blocks(:,2)==1);

for ii=1:size(yr_change,1)
    if ii<size(yr_change,1)
        clim_blocks(yr_change(ii):yr_change(ii+1)-1,3)=ii+climate_20thCent(end,3);
    else
        clim_blocks(yr_change(ii):end,3)=ii+climate_20thCent(end,3);
    end
% clim_blocks(:,3)=repmat(climate_20thCent(:,3),2,1)
climate=[climate_spinup;climate_20thCent;clim_blocks];
end
end
