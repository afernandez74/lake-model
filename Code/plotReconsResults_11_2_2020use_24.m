clear

%%
plotPT20thCent_24
dates_P_seasonality = dates_years;
%%
%load level logger lake level
LLobs=readmatrix('../Data/castorLakeLevel0521.csv');
dates_obs=datetime(LLobs(:,3),LLobs(:,2),LLobs(:,1));
%load satellite image lake level
LLsat=readmatrix('../Data/SatDepths.csv');
dates_sat=datetime(LLsat(:,3),LLsat(:,2),LLsat(:,1));
%load 18O sample values
d18O_samples=readmatrix('../Data/18O_samples_Castor_2.csv');
dates_samples=datetime(d18O_samples(:,3),d18O_samples(:,2),d18O_samples(:,1));
%load Sam's annual 100-year 18O record
d18O_core=readmatrix('../Data/18OannualSam.csv');
dates_core=datetime(d18O_core(:,1),07,01);

% load SUG tree ring width data
load('../Data/SUG_TRWdata_20thCent.mat');
dates_SUG=datetime(SUG_TRWdata_20thCent(:,1),07,01);

%load monte carlo simulation results 
load('../Results/MCres_recons_ClimMDaragTiming_21October_2024_12 42_i=1000.mat')
results1=results;
load('../Results/MCres_recons_ClimMDaragTiming(narrowClimRanges_seasonalTemp)_21October_2024_16 04_i=1000.mat')
results2=results;
clear 'results'
%%

daily_LD_all=[results1.daily_LD_mm results2.daily_LD_mm];
[~,cols_w_nan]=find(isnan(daily_LD_all));
cols_w_nan=unique(cols_w_nan);
daily_LD_all(:,cols_w_nan)=[];
mean_daily_LD=mean(daily_LD_all,2);
std_daily_LD=std(daily_LD_all,0,2);

figure

plot(results1.dates,mean_daily_LD,'k-','LineWidth',1);
hold
plot(results1.dates,mean_daily_LD+2*std_daily_LD,'k--','LineWidth',1)
plot(results1.dates,mean_daily_LD-2*std_daily_LD,'k--','LineWidth',1)
plot(dates_obs,LLobs(:,4),'r-','LineWidth',0.8)
plot(dates_sat,LLsat(:,4),'b.','Markersize',20)



ylabel('Lake Depth (cm)', 'FontSize', 20)
grid on
title('20th Century Modeled and Observed Lake Depth', 'FontSize', 25)
ax=gca;
ax.FontSize=20;
% ax.XLim=[1 12];
% ax.XTickLabel={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ax.XMinorTick = 'off';


legend({'Modeled Daily Water Depth','Measured Lake Depth (levelogger)','Estimated Lake Level (satelite images)'},...
    'FontSize', 15,'Location','southeast')
%%
% compute 10 year trailing averages for the yearly arag 18O records

goal=-2.4;
uncertainty=0.2;

dates_years=results1.dates_years;
dates_years_20thCent=dates_years(20:find(dates_years==datetime(2021,1,1)));
dates_years_no_mod=dates_years(find(dates_years==datetime(2021,1,1))+20:find(dates_years==datetime(2021,1,1))+118);
dates_years_mod=dates_years(end-100:end);

start_20thCent=find(year(results1.dates_years)==1900);
end_20thCent=find(year(results1.dates_years)==2019);

year_arag_all=[results1.trail_avg_year_summer_arag results2.trail_avg_year_summer_arag];
mean_arag_20thCent_all=[results1.mean_arag_20thCent results2.mean_arag_20thCent];
mean_arag_no_mod_all=[results1.mean_arag_no_mod results2.mean_arag_no_mod];
mean_arag_mod_all=[results1.mean_arag_mod results2.mean_arag_mod];

[~,cols_w_nan]=find(isnan(year_arag_all));
cols_w_nan=unique(cols_w_nan);
year_arag_all(:,cols_w_nan)=[];
mean_arag_20thCent_all(:,cols_w_nan)=[];
mean_arag_no_mod_all(:,cols_w_nan)=[];
mean_arag_mod_all(:,cols_w_nan)=[];

mean_SamCore=mean(d18O_core(:,2));
pos_good=find(mean_arag_20thCent_all < mean_SamCore+uncertainty & mean_arag_20thCent_all > mean_SamCore-uncertainty);
pos_good_recons=find(mean_arag_mod_all < goal+uncertainty*2 & mean_arag_mod_all > goal-uncertainty*2);
mean_year_arag_all=mean(year_arag_all(:,pos_good),2);
std_year_arag_all=std(year_arag_all(:,:),0,2);
mean_20thCentModel=mean(mean_year_arag_all(start_20thCent:end_20thCent));

year_arag_all_20thCent=year_arag_all(20:find(dates_years==datetime(2021,1,1)),:);
year_arag_all_no_mod=year_arag_all(find(dates_years==datetime(2021,1,1))+20:find(dates_years==datetime(2021,1,1))+118,:);
year_arag_all_mod=year_arag_all(end-100:end,:);

figure
hold

% for i=1:size(year_arag_all,2)
%     plot(results1.dates_years,year_arag_all(:,i),'-','Color',[.9 .9 .9])
% end
mean_hist = mean(year_arag_all_20thCent(:,pos_good_recons),2);
std_hist = std(year_arag_all_20thCent(:,pos_good_recons),0,2);

%plot(dates_years_20thCent,year_arag_all_20thCent(:,pos_good_recons),'-','Color',[.9 .9 .9])
plot(dates_years_20thCent,mean_hist,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_years_20thCent,mean_hist+std_hist,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_years_20thCent,mean_hist-std_hist,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_core,movmean(d18O_core(:,2), [10 0]),'k-.','MarkerSize',10,'LineWidth',3)

plot([results1.dates_years(start_20thCent),results1.dates_years(end_20thCent)],[mean_20thCentModel mean_20thCentModel],':','LineWidth',3,'Color',[.4 .4 .4])
plot([dates_core(1),dates_core(end)],[mean_SamCore mean_SamCore],'k:','LineWidth',3)

legend({'1000 realization mean 10-yr trailing avg modeled sed \delta ^{18}O','Measured Sediment Core \delta ^{18}O','Model 20th Century mean value ','Core 20th Century mean value'},...
    'FontSize', 15,'Location','southwest')

ax=gca;
ax.YLabel.String=(strcat('VPDB \delta ^{18}O ','( ',char(8240),')'));
ax.XMinorTick = 'off';
ax.XLabel.String=('Year');
ax.FontSize=20;
title ('Castor Lake Modeled and Observed Sediment \delta^{18}O over the 20th Century','FontSize',30)

%% plot same figure as above but with the no_mod and mod stuffs
figure
hold
mean_no_mod = mean(year_arag_all_no_mod(:,pos_good_recons),2);
std_no_mod = std(year_arag_all_no_mod(:,pos_good_recons),0,2);

plot(dates_years_no_mod,year_arag_all_no_mod(:,pos_good_recons),'-','Color',[.9 .9 .9])
plot(dates_years_no_mod,mean_no_mod+std_no_mod,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_years_no_mod,mean_no_mod-std_no_mod,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_years_no_mod,mean_no_mod,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot([dates_years_no_mod(1),dates_years_no_mod(end)],[mean(mean(year_arag_all_no_mod(:,pos_good_recons),2)) mean(mean(year_arag_all_no_mod(:,pos_good_recons),2))],'k:','LineWidth',3)

ax=gca;
ax.YLabel.String=(strcat('VPDB \delta ^{18}O ','( ',char(8240),')'));
ax.XMinorTick = 'off';
ax.XLabel.String=('Year');
ax.FontSize=20;

figure
hold

mean_mod = mean(year_arag_all_mod(:,pos_good_recons),2);
std_mod = std(year_arag_all_mod(:,pos_good_recons),0,2);

plot(dates_years_mod,year_arag_all_mod(:,pos_good_recons),'-','Color',[.9 .9 .9])
plot(dates_years_mod,mean_mod,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_years_mod,mean_mod + std_mod,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot(dates_years_mod,mean_mod - std_mod,'k.-','LineWidth',2,'Color',[.4 .4 .4]);
plot([dates_years_mod(1),dates_years_mod(end)],[mean(mean(year_arag_all_mod(:,pos_good_recons),2)) mean(mean(year_arag_all_mod(:,pos_good_recons),2))],'k:','LineWidth',3)

ax=gca;
ax.YLabel.String=(strcat('VPDB \delta ^{18}O ','( ',char(8240),')'));
ax.XMinorTick = 'off';
ax.XLabel.String=('Year');
ax.FontSize=20;


%%
%plot same data as above, with SUG TRD data all as zcores
figure
mean_year_arag_all_20thCent=mean(year_arag_all_20thCent(:,pos_good_recons),2);
plot(dates_years_20thCent,zscore(mean_year_arag_all_20thCent),'k-','LineWidth',1);
hold
plot(dates_core,zscore(movmean(d18O_core(:,2),[10 0])),'b-','MarkerSize',10,'LineWidth',0.5)
plot(dates_SUG,-1*zscore(movmean(SUG_TRWdata_20thCent(:,2),[10 0])),'r-','MarkerSize',10,'LineWidth',0.5)
plot(dates_P_seasonality,-1*zscore(movmean(year_winter_p, [10 0])),'-','linewidth',0.5,'color',[.5 .5 .5])
% go to plotPT20thCent for winter precip data 
legend({'1000 realization mean 10-yr trailing avg modeled sed \delta ^{18}O','Measured Sediment Core \delta ^{18}O',...
    'Total Ring Width SUG lake 10-yr moving avg','winter precip'},...
    'FontSize', 14,'Location','northwest')

ax=gca;
ax.YLabel.String=('zscore');
ax.XMinorTick = 'on';
ax.XLabel.String=('Year');
ax.FontSize=20;
grid on 
title ({'Castor Lake Modeled and Observed Sediment \delta ^{18} O';'and TRW from SUG over the 20th Century'},'FontSize',30)

%%
% figure 
% hold
% 
% for i=1:size(results1.MD_max_v,2)
%     plot(results1.dates,results1.daily_dl(:,i),'LineWidth',0.1);
% end
% 
% plot(dates_samples,d18O_samples(:,4),'o','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
% ax=gca;
% ax.FontSize=20;
% ax.XMinorTick = 'off';