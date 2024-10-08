%% monte carlo simulations
clear
tic
% created from 20th century uncertainty calibration in
% calibrationResults20thCentUncertCalib.m
load('../Data/uncertainty_final.mat');

clim=create20thCenturyClimate_Spinup_24(0);
final_20thCentYear=clim(end,3);
years_clim=(clim(1,3):clim(end,3))';
n_spinupyears=20; %number of years average climate spinup of the model
n_reps=2; %number of times the 20th Cent dataset will be permutated in its re-ordered version by 10yr blocks

climate=create20thCenturyClimate_Spinup_20thCent_RandYrBlocks_24(n_spinupyears,n_reps);
    
results=struct;
ni_MC=1000;
a=2.45;
b=38.03;
C_IN=0.001;
AWC_mod=1.87;

P_mod_summer_mu=-0.0144;
P_mod_summer_sigma=0.19*2;

P_mod_winter_mu=-0.1874;
P_mod_winter_sigma=0.0574*2;

T_mod_summer_mu=0.0;
T_mod_summer_sigma=0.5*2;
T_mod_summer_range=[-2 -2];

T_mod_winter_mu=-0.0;
T_mod_winter_sigma=0.5*2;
T_mod_winter_range=[-2 2];

RH_mod_mu=-0.188;
RH_mod_sigma=0.0531*2;

results.MD_max_v=NaN(1,ni_MC);
results.summer_begin_v=NaN(1,ni_MC);
results.summer_len_v=NaN(1,ni_MC);
results.P_mod_summer_v=NaN(1,ni_MC);
results.P_mod_winter_v=NaN(1,ni_MC);
results.T_mod_summer_v=NaN(1,ni_MC);
results.T_mod_winter_v=NaN(1,ni_MC);
results.RH_mod_v=NaN(1,ni_MC);

results.year_summer_arag=NaN(size((climate(1,3):climate(end,3))',1),ni_MC);
results.trail_avg_year_summer_arag=NaN(size((climate(1,3):climate(end,3))',1),ni_MC);

results.obs_coefficients=NaN(1,ni_MC); %intercept row 1, %slope row 2
results.obs_corr=NaN(1,ni_MC);
results.obs_rmse=NaN(1,ni_MC);
results.obs_mean_diff=NaN(1,ni_MC);

results.samples_rmse=NaN(1,ni_MC);

results.core_corr=NaN(1,ni_MC);
results.core_rmse=NaN(1,ni_MC);
results.core_rmse_mm=NaN(1,ni_MC);

results.meanCore_meanYearArag_diff=NaN(1,ni_MC);
% % results=zeros(ni_MC,5);

for i_MC = 1 : ni_MC
        
    MD_max=MD_new_mean+randn*MD_new_std;
    if MD_max <= 10
        MD_max=10;
    end
    results.MD_max_v(i_MC)=MD_max;
    
    summer_begin=floor(summer_begin_new_mean+randn*summer_begin_new_std);
    if summer_begin < 105
        summer_begin=105;
    end
    if summer_begin > 181
        summer_begin=181;
    end
    results.summer_begin_v(i_MC)=summer_begin;
    
    summer_len=floor(summer_len_new_mean+randn*summer_len_new_std);
    results.summer_len_v(i_MC)=summer_len;
    
    P_mod_summer=(P_mod_summer_mu+randn*P_mod_summer_sigma);
    results.P_mod_summer_v(i_MC)=P_mod_summer;
    
    P_mod_winter=(P_mod_winter_mu+randn*P_mod_winter_sigma);
    results.P_mod_winter_v(i_MC)=P_mod_winter;
    
    T_mod_summer=(T_mod_summer_mu+randn*T_mod_summer_sigma);
    results.T_mod_summer_v(i_MC)=T_mod_summer;
    
    T_mod_winter=(T_mod_winter_mu+randn*T_mod_winter_sigma);
    results.T_mod_winter_v(i_MC)=T_mod_winter;
    
    RH_mod=(RH_mod_mu+randn*RH_mod_sigma);
    results.RH_mod_v(i_MC)=RH_mod;
    
    yearly_P_seasonality=calc_yearly_P_seasonality_indep(P_mod_summer,P_mod_winter,climate);
    ModelCastorContClimDaily_20thCenturyClimMC_ReconstLoops_24;
    
    results.daily_LD_mm(:,i_MC)=daily_LDv_mm;
    results.daily_dl(:,i_MC)=daily_dl;
    
    if size(year_summer_arag,1)==size(results.trail_avg_year_summer_arag,1)
        results.year_summer_arag(:,i_MC)=year_summer_arag;
        results.trail_avg_year_summer_arag(:,i_MC)=trail_avg_year_summer_arag;
    else 
        results.trail_avg_year_summer_arag(:,i_MC)=NaN;
    end
    
    
    if bad_run==false
        
        results.obs_coefficients(1:2,i_MC)=lm_obs.Coefficients.Estimate;
        results.obs_corr(i_MC)=obs_corr;
        results.obs_rmse(i_MC)=obs_rmse;
        results.obs_mean_diff(i_MC)=obs_mean_diff;

        results.samples_rmse(i_MC)=obs_rmse;

        results.core_corr(i_MC)=core_corr;
        results.core_rmse(i_MC)=core_rmse;
        results.meanCore_meanYearArag_diff(i_MC)=meanCore_meanYearArag_diff;
    end
    
    i_MC

    climate=create20thCenturyClimate_Spinup_20thCent_RandYrBlocks_24(n_spinupyears,n_reps);
end

total_years=size(results.trail_avg_year_summer_arag,1);
total_years_v=(1:total_years)';
results.dates=dates_clim;
results.dates_years=dates_years;
yearsi=year(dates_years);

results.mean_arag_spinup=mean(results.year_summer_arag(1:n_spinupyears,:),1);
% results.mean_arag_20thCent=mean(results.year_arag(n_spinupyears+20:end,:),1);
results.mean_arag_20thCent=mean(results.year_summer_arag(n_spinupyears+20:find(yearsi==2018),:),1);
results.mean_arag_no_mod=mean(results.year_summer_arag(find(yearsi==2018)+20:find(yearsi==floor(mean([2018 2256]))),:),1);
results.mean_arag_mod=mean(results.year_summer_arag(find(yearsi==floor(mean([2018 2256])))+20:end,:),1);


%% stuff for naming file


filename=strcat('../Results/MCres_recons_ClimMDaragTiming(narrowClimRanges_seasonalTemp)'...
    ,datestr(now,'_ddmmmm_yyyy_HH MM'),'_i=',num2str(ni_MC),'.mat');

save(filename,'results','-v7.3');
toc

%% reconstruction

goal=-2.4;
uncertainty=0.2;
positions=find(results.mean_arag_mod < goal+uncertainty & results.mean_arag_mod > goal-uncertainty);
arag_recons=results.mean_arag_mod(positions);
MD_recons=results.MD_max_v(positions);
summer_begin_recons=results.summer_begin_v(positions);
summer_len_recons=results.summer_len_v(positions);
P_mod_summer_recons=results.P_mod_summer_v(positions);
P_mod_winter_recons=results.P_mod_winter_v(positions);
T_mod_summer_recons=results.T_mod_summer_v(positions);
T_mod_winter_recons=results.T_mod_winter_v(positions);
RH_mod_recons=results.RH_mod_v(positions);

figure
subplot(4,1,1)
histogram(results.mean_arag_mod,[-8:0.2:0])
title('All arag results')

subplot(4,1,2)
histogram(MD_recons)
title('MD reconstruction values')

subplot(4,1,3)
histogram(summer_begin_recons)
title('beginning of summer')


subplot(4,1,4)
histogram(summer_len_recons)
title('Length of summer')
%% 
load('D:\MCres_Castor\20thCentReconstResults(loopsClim&Uncert)\V_3\MCres_recons_ClimMDaragTiming(narrowClimRanges_seasonalTemp)_02November_2020_21 35_i=1000.mat')

mean_P_mod_summer_recons=mean(P_mod_summer_recons);
mean_P_mod_winter_recons=mean(P_mod_winter_recons);
mean_T_mod_summer_recons=mean(T_mod_winter_recons);
mean_T_mod_winter_recons=mean(T_mod_winter_recons);
mean_RH_mod_recons=mean(RH_mod_recons);
sigma_P_mod_summer_recons=[mean_P_mod_summer_recons-std(P_mod_summer_recons)...
     mean_P_mod_summer_recons+std(P_mod_summer_recons)];
sigma_P_mod_winter_recons=[mean_P_mod_winter_recons-std(P_mod_winter_recons)...
     mean_P_mod_winter_recons+std(P_mod_winter_recons)];
sigma_T_mod_summer_recons=[mean_T_mod_summer_recons-std(mean_T_mod_summer_recons)...
     mean_T_mod_summer_recons+std(mean_T_mod_summer_recons)];
sigma_T_mod_winter_recons=[mean_T_mod_winter_recons-std(mean_T_mod_winter_recons)...
     mean_T_mod_winter_recons+std(mean_T_mod_winter_recons)];
sigma_RH_mod_recons=[mean_RH_mod_recons-std(RH_mod_recons)...
     mean_RH_mod_recons+std(RH_mod_recons)];

figure;hold;

subplot(5,2,1)
ksdensity(P_mod_summer_recons)
xline(mean(P_mod_summer_recons))
xline(mean(P_mod_summer_recons)+std(P_mod_summer_recons))
xline(mean(P_mod_summer_recons)-std(P_mod_summer_recons))
xlim([-0.5 0.5])
title('Summer precip modification')

subplot(5,2,2)
cdfplot(P_mod_summer_recons)
xlim([-0.5 0.5])
title('Summer precip modification ECDF')

subplot(5,2,3)
ksdensity(P_mod_winter_recons)
xline(mean(P_mod_winter_recons))
xline(mean(P_mod_winter_recons)+std(P_mod_winter_recons))
xline(mean(P_mod_winter_recons)-std(P_mod_winter_recons))
xlim([-0.5 0.5])
title('Winter precip modification')

subplot(5,2,4)
cdfplot(P_mod_winter_recons)
xlim([-0.5 0.5])
title('Winter precip modification ECDF')

subplot(5,2,5)
ksdensity(T_mod_winter_recons)
xline(mean(T_mod_winter_recons))
xline(mean(T_mod_winter_recons)+std(T_mod_winter_recons))
xline(mean(T_mod_winter_recons)-std(T_mod_winter_recons))
xlim([-4 4])
title('Temp winter modification (^{o}C)')

subplot(5,2,6)
cdfplot(T_mod_winter_recons)
xlim([-4 4])
title('Temp winter modification (^{o}C) ECDF')

subplot(5,2,7)
ksdensity(T_mod_summer_recons)
xline(mean(T_mod_summer_recons))
xline(mean(T_mod_summer_recons)+std(T_mod_summer_recons))
xline(mean(T_mod_summer_recons)-std(T_mod_summer_recons))
xlim([-4 4])
title('Temp summer modification (^{o}C)')

subplot(5,2,8)
cdfplot(T_mod_summer_recons)
xlim([-4 4])
title('Temp summer modification (^{o}C) ECDF')

subplot(5,2,9)
ksdensity(RH_mod_recons)
xline(mean(RH_mod_recons))
xline(mean(RH_mod_recons)+std(RH_mod_recons))
xline(mean(RH_mod_recons)-std(RH_mod_recons))
xlim([-0.5 0.5])
title('RH modification')

subplot(5,2,10)
cdfplot(RH_mod_recons)
xlim([-0.5 0.5])
title('RH modification')

%% all centered climate 
load('D:\MCres_Castor\20thCentReconstResults(loopsClim&Uncert)\V_3\MCres_recons_ClimMDaragTiming_31October_2020_08 09_i=1000.mat')
P_mod_summer_recons=results.P_mod_summer_v;
P_mod_winter_recons=results.P_mod_winter_v;
RH_mod_recons=results.RH_mod_v;
figure;hold;
T_mod_summer_range=[-2 2];
T_mod_winter_range=[-2 2];

T_mod_summer_mu=mean(T_mod_summer_range);
T_mod_summer_sigma=(T_mod_summer_range(2)-T_mod_summer_mu)/2;

T_mod_winter_mu=mean(T_mod_winter_range);
T_mod_winter_sigma=(T_mod_winter_range(2)-T_mod_winter_mu)/2;

T_mod_summer_recons=T_mod_summer_mu+randn(1000,1)*T_mod_summer_sigma;

T_mod_winter_recons=T_mod_winter_mu+randn(1000,1)*T_mod_winter_sigma;


subplot(5,2,1)
ksdensity(P_mod_summer_recons)
xline(mean(P_mod_summer_recons))
xline(mean(P_mod_summer_recons)+std(P_mod_summer_recons))
xline(mean(P_mod_summer_recons)-std(P_mod_summer_recons))
xlim([-0.5 0.5])
title('Summer precip modification')

subplot(5,2,3)
ksdensity(P_mod_winter_recons)
xline(mean(P_mod_winter_recons))
xline(mean(P_mod_winter_recons)+std(P_mod_winter_recons))
xline(mean(P_mod_winter_recons)-std(P_mod_winter_recons))
xlim([-0.5 0.5])
title('Winter precip modification')

subplot(5,2,5)
ksdensity(T_mod_winter_recons)
xline(mean(T_mod_winter_recons))
xline(mean(T_mod_winter_recons)+std(T_mod_winter_recons))
xline(mean(T_mod_winter_recons)-std(T_mod_winter_recons))
xlim([-4 4])
title('Temp winter modification (^{o}C)')

subplot(5,2,7)
ksdensity(T_mod_summer_recons)
xline(mean(T_mod_summer_recons))
xline(mean(T_mod_summer_recons)+std(T_mod_summer_recons))
xline(mean(T_mod_summer_recons)-std(T_mod_summer_recons))
xlim([-4 4])
title('Temp summer modification (^{o}C)')

subplot(5,2,9)
ksdensity(RH_mod_recons)
xline(mean(RH_mod_recons))
xline(mean(RH_mod_recons)+std(RH_mod_recons))
xline(mean(RH_mod_recons)-std(RH_mod_recons))
xlim([-0.5 0.5])
title('RH modification')
