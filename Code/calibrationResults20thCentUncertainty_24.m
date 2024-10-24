%%
clear
load ('../Results/MCres_recons_MDaragTiming_18October_2024_16 53_i=1000first.mat')
%load('D:\MCres_Castor\20thCentUncertCalib(MD and summer timing)\V_3\MCres_recons_MDaragTiming_30October_2020_16 32_i=1000first.mat')
%%
%load 18O sample values
d18O_samples=readmatrix('../Data/18O_samples_Castor_2.csv');
dates_samples=datetime(d18O_samples(:,3),d18O_samples(:,2),d18O_samples(:,1));
%load Sam's annual 100-year 18O record
d18O_core=readmatrix('../Data/18OannualSam.csv');
dates_core=datetime(d18O_core(:,1),07,01);
%%
climate=create20thCenturyClimate_Spinup_24(20);
years=(climate(1,3):climate(end,3))';
dates_clim=datetime(climate(:,3),climate(:,2),climate(:,1));
dates_clim_years=datetime(years(:),07,01);
dates_clim_years=unique(dates_clim_years);
dates_clim_years(isnat(dates_clim_years))=[];
%%
mean_core=mean(d18O_core(:,2));
std_core=std(d18O_core(:,2));

range_20thCent=[mean_core - 0.2   mean_core + 0.2];
% year_arag_trailavg=movmean(results.year_arag,[10 0],1);
mean_year_arag_trailavg=mean(results.trail_avg_year_summer_arag,1);

pos_good=mean_year_arag_trailavg>=range_20thCent(1) & mean_year_arag_trailavg<=range_20thCent(2);

MD_good=results.MD_max_v(pos_good);
MD_new_mean=mean(MD_good);
MD_new_std=std(MD_good);

summer_begin_good=results.summer_begin_v(pos_good);
summer_begin_new_mean=mean(summer_begin_good);
summer_begin_new_std=std(summer_begin_good);

summer_len_good=results.summer_len_v(pos_good);
summer_len_new_mean=mean(summer_len_good);
summer_len_new_std=std(summer_len_good);

%%

save('../Data/uncertainty_final_24','MD_new_mean','MD_new_std','summer_begin_new_mean','summer_begin_new_std','summer_len_new_mean','summer_len_new_std');