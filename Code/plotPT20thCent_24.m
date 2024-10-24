clear

%load ('G:\My Drive\lake_model\LakeModel2\FinalClimateCompilation\climateCompilation20thCent.mat')
load('../Data/climateCompilation20thCent_24.mat')

%% calculate yearly total precipitation amounts and average temperatures
year_precip=NaN(final_pt(end,3)-final_pt(1,3)-1,1);
year_temp=NaN(final_pt(end,3)-final_pt(1,3)-1,1);
j=1;

for i = final_pt(1,3)+1:final_pt(end,3)-1
    
    year_precip(j)=sum(final_pt(final_pt(:,3)==i,4));
    year_temp(j)=mean(final_pt(final_pt(:,3)==i,5));
    j=j+1;
    
end
dates_years=(final_pt(1,3)+1):(final_pt(end,3)-1);
dates_years=(datetime(dates_years,07,01))';

%% calculate warm season vs cold season yearly precipitation amounts 

% winter starts end of sept
towinter=[9,30];
%summer starts end of march
tosummer=[3,31];

temp_winter_p=0.0;
temp_summer_p=0.0;
j=1;
year_summer_p=NaN(final_pt(end,3)-final_pt(1,3)-1,1);
year_winter_p=NaN(final_pt(end,3)-final_pt(1,3)-1,1);

for i=1:size(final_pt,1)
    if final_pt(i,2)<4 || final_pt(i,2)>=10
        temp_winter_p=temp_winter_p+final_pt(i,4);
    else
        temp_summer_p=temp_summer_p+final_pt(i,4);
    end
    if [final_pt(i,2) final_pt(i,1)] == towinter
        year_summer_p(j)=temp_summer_p;
        temp_summer_p=0;
        j=j+1;
    elseif [final_pt(i,2),final_pt(i,1)] == tosummer
        year_winter_p(j)=temp_winter_p;
        temp_winter_p=0;
    end
    if i==size(final_pt,1)
        year_summer_p(end)=temp_summer_p;
    end        
end
year_winter_p(1)=[];



%% use tight plot

ha = tight_subplot(2,1,[.005 .5],[.05 .05],[.1 .1]);
axes(ha(1));
hold on 
plot(dates_years,movmean(year_precip,[2.5 2.5]),'-ok','linewidth',1)
% plot(dates_years,year_summer_p+year_winter_p,'--k','LineWidth',0.5)
plot(dates_years,movmean(year_summer_p,[2.5 2.5]),'--s','Color',[0.25 0.25 0.25],'LineWidth',0.1)
plot(dates_years,movmean(year_winter_p,[2.5 2.5]),'--*','Color',[0.5 0.5 0.5],'LineWidth',0.1)
grid on 
ax=gca;
ax.FontSize=23;
ax.YLabel.String='Annual Precipitation Amount (mm)';
ax.XTickLabel=[];
ax.YTickLabelMode='auto';
yticks([100 300 500])
box on
legend ('Total Year Precipitation','Year Warm Season Precipitation','Year Cold Season Precipitation')
ylim([0 600])
xlim([dates_years(1) dates_years(end)])

axes(ha(2));
plot(dates_final_pt,movmean(final_pt(:,5),360),'-','LineWidth',0.1,'Color',[0.4 0.4 0.4])
hold
plot(dates_years,movmean(year_temp,[2 2]),'-k','LineWidth',1)
grid on
ax=gca;
ax.FontSize=23;
ax.YLabel.String='Average Yearly Air Temperature';
ax.YAxisLocation='right';
yticks([4 6 8 10])
ylim([4 11])
xlim([dates_years(1) dates_years(end)])

ylh = get(gca,'ylabel');
set(ylh,'rotation',270,'VerticalAlignment','middle')
legend('Yearly Average Temperature','5 year moving mean Temperature')
