%% monte carlo simulations
clear
tic

climate=create20thCenturyClimate_Spinup_24(20);

results=struct;
ni_MC=1000;
a=2.45;
b=38.03;
C_IN=0.001;
AWC_mod=1.87;
MD_max_mu=600;
MD_max_sigma=MD_max_mu/4;
% MD_max_sigma=0;

summer_begin_range=[105 181];
summer_begin_mu=floor(mean(summer_begin_range));
summer_begin_sigma=floor((summer_begin_range(2)-summer_begin_mu)/2);
% summer_begin_sigma=0;

summer_len_range=[77 122];
summer_len_mu=floor(mean(summer_len_range));
summer_len_sigma=floor((summer_len_range(2)-summer_len_mu)/2);
% summer_len_sigma=0;

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

for i_MC = 1 : ni_MC
    MD_max=MD_max_mu+randn*MD_max_sigma;
    if MD_max <= 10
        MD_max=10;
    end
    results.MD_max_v(i_MC)=MD_max;
    
    summer_begin=floor(summer_begin_mu+randn*summer_begin_sigma);
    if summer_begin < 105
        summer_begin=105;
    end
    if summer_begin > 181
        summer_begin=181;
    end
    results.summer_begin_v(i_MC)=summer_begin;
    
    summer_len=floor(summer_len_mu+randn*summer_len_sigma);
    results.summer_len_v(i_MC)=summer_len;
    
    ModelCastorContClimDaily_20thCenturyClimateMC_UncertCalib_24;
    
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
end


%% stuff for naming file


% filename=strcat('D:\MCres_Castor\20thCentUncertCalib(MD and summer timing)\V_3\MCres_recons_MDaragTiming'...
%     ,datestr(now,'_ddmmmm_yyyy_HH MM'),'_i=',num2str(ni_MC),'first.mat');

filename = strcat('../Results/MCres_recons_MDaragTiming'...
     ,datestr(now,'_ddmmmm_yyyy_HH MM'),'_i=',num2str(ni_MC),'first.mat');


save(filename,'results','-v7.3');
toc

% filename=convertCharsToStrings(filename);

save(filename,'results');
toc