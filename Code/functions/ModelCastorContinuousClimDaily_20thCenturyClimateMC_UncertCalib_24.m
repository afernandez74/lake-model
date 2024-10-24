%Initialize climate_20thCent parameters
%run createContinuousDailyDataXXX.m first%

climate=create20thCenturyClimate_Spinup_24(20);

Ta_year=createDailyInterpTdataFromAVG();
WS_year=createDailyInterpWSdataFromAVG;
Ra_year=createDailyInterpRadataFromAVG;
P_d18O_year=createDailyInterpP_d18OdataFromAVG;
createMD_mod;
MD_year=dailyAVGMD;

% Set simulation parameters 
t0=0;                           % Initial time [days]
tf=size(climate,1);             % Final time [days]
dt=1/(2^1);                     % Time step [days]
ni=(tf-t0)/dt+1;                % Number of iterations []

MIX_OFF=false;


% Lake catchment constants

AWC=4.6/100;                    % total available water capacity
AWC=AWC*AWC_mod;                % soil moisture modifying factor for
                                % volcanic effects
AWC_SSf=0.5;                    % Proportion surface soil / deep soil
AWC_SS=AWC*AWC_SSf;             % Available water capacity in the surface soil layer [m]
AWC_DS=AWC*(1-AWC_SSf);         % Available water capacity in the deep soil layer [m]

TCA=0.86 * 1e+06;               % Total catchment area [m2]
ALB_L=0.08;                     % Albedo Lake []
ALB_G=0.25;                     % Albedo Catchment []
               

% C_IN=0.0035;
C = 14.3;                       % Experimentally derived isotopic separation value
Outflow_delay=0.05;

% Loads Hypsometry data from Hyps function
% constants for exponential outseepage function

load 'CastorNewBathymetry.mat';
HypsOutseep = calcOutSeepageData_daily(a,b,castorBathymetry);
Max_depth=1335;
Max_V=interp1(HypsOutseep(:,1),HypsOutseep(:,3),Max_depth); 

SLi=216142;
DLi=77445;
SSi=0;
DSi=0;
INi=50000;
SPi=0;
LVi=SLi+DLi;

SL_isoi=-172914;
DL_isoi=-61956;
SS_isoi=0;
DS_isoi=0;
IN_isoi=-150000;
SP_isoi=0;

% Set initial conditions for the six hydro reservoirs
SL=NaN(ni,1);
DL=NaN(ni,1);
SS=NaN(ni,1);
DS=NaN(ni,1);
IN=NaN(ni,1);
SP=NaN(ni,1);

SL(1)=SLi;
DL(1)=DLi;
SS(1)=SSi;
DS(1)=DSi;
IN(1)=INi;
SP(1)=SPi;
% Set initial conditions for the isotope reservoirs

SL_iso=NaN(ni,1);
DL_iso=NaN(ni,1);
SS_iso=NaN(ni,1);
DS_iso=NaN(ni,1);
IN_iso=NaN(ni,1);
SP_iso=NaN(ni,1);

SL_iso(1)=SL_isoi;
DL_iso(1)=DL_isoi;
SS_iso(1)=SS_isoi;
DS_iso(1)=DS_isoi;
IN_iso(1)=IN_isoi;
SP_iso(1)=SP_isoi;

% Allocate space for plotting different fluxes and other variables

T=NaN(ni,1);
Vv=NaN(ni,1);
Vv(1)=LVi;
SAv=NaN(ni,1);
LDv=NaN(ni,1);
CAv=NaN(ni,1);
fpv=NaN(ni,1);
finv=NaN(ni,1);
fdlmv=NaN(ni,1);
fslmv=NaN(ni,1);
frov=NaN(ni,1);
fdsdv=NaN(ni,1);
frv=NaN(ni,1);
fdosv=NaN(ni,1);
fsosv=NaN(ni,1);
fdsev=NaN(ni,1);
fev=NaN(ni,1);
fsfv=NaN(ni,1);
fsmv=NaN(ni,1);
fssdv=NaN(ni,1);
fssev=NaN(ni,1);
fssiv=NaN(ni,1);
de=NaN(ni,1);
dl=NaN(ni,1);
svcv=NaN(ni,1);
fofv=NaN(ni,1);
frssv=NaN(ni,1);
frinv=NaN(ni,1);
fsmssv=NaN(ni,1);
fsminv=NaN(ni,1);
petdemv=NaN(ni,1);

years=(climate(1,3):climate(end,3))';
dates_years=datetime(years,1,1);
year_winter_p=NaN(size(years,1)-1,1);
year_summer_p=NaN(size(years,1)-1,1);

summer_arag=NaN(ni,1);
summer_dl=NaN(ni,1);
d18O_aragonite=NaN(ni,1);
year_summer_arag=NaN(ni,1);
year_summer_dl=NaN(ni,1);
dayi=NaN;monthi=NaN;year=NaN;
totaldays=(1:size(climate,1))';
i=1;
Tw_switch=0;
loop_counter=-1;
summer_end=summer_begin+summer_len;
if summer_end > 227
    summer_end=227;
    summer_len=summer_end-summer_begin;
end
   
midsummer=floor((summer_begin+summer_end)/2);
temp_summer_arag=NaN(summer_len/dt,1);
temp_summer_dl=NaN(summer_len/dt,1);
days=NaN(ni,1);
months=NaN(ni,1);
years=NaN(ni,1);

i_summer=0;

sizedat=size(climate,1);

for t=t0:dt:tf-1

%apply mod sizedat floor t so it loops after the 12 years
T(i)=t;
if mod((t),sizedat)==0
    loop_counter=loop_counter+1;
end

dayi=climate(mod(floor(t),sizedat)+1,1);
days(i)=dayi;
monthi=climate(mod(floor(t),sizedat)+1,2);
months(i)=monthi;
yearinloop=climate(mod(floor(t),sizedat)+1,3);
yeari=climate(mod(floor(t),sizedat)+1,3)+12*loop_counter;
years(i)=yeari;

dayinyear=day365(monthi,dayi);

Pi = climate(mod(floor(t),sizedat)+1,4);        %[mm/day]
Tai = climate(mod(floor(t),sizedat)+1,5);       %[C]
RHi = climate(mod(floor(t),sizedat)+1,6);       %[%]
Rsi = climate(mod(floor(t),sizedat)+1,7);       %[MJ/day/m2] insolation
Twi = climate(mod(floor(t),sizedat)+1,8);
MDi = MD_year(dayinyear);                                %[cm]
Rai = Ra_year(dayinyear);                                %[MJ/day/m2] extraterrestrial solar rad
WSi = WS_year(dayinyear);                                 %[m/s]

Tai_Tavg_diff=Tai-Ta_year(dayinyear);
P_d18Oi = P_d18O_year(dayinyear)+0.6*Tai_Tavg_diff;  %[permil]

SAi = SA(LVi,castorBathymetry);                 % Lake Surface Area at iteration [m2]
LDi = LD(LVi,castorBathymetry);LDv(i)=LDi;      % Lake Depth at iteration [cm]


CAi = TCA - SAi;                    % Catchment Area at iteration [m2]
DLDi = LD(LVi,castorBathymetry) - MDi;          % Deep Lake Depth at iteration [cm]
DLVi = DLV(DLDi,castorBathymetry);              % Deep Lake Volume at iteration [m3] 
SVCi = LVi - DLVi;                  % Surface Volume Control at iteration [m3]
svcv(i)=SVCi;
SSC = CAi * AWC_SS;                 % Surface Soil Capacity [m3]
DSC = CAi * AWC_DS;                 % Deep Soil Capacity [m3]
TSC = SSC + DSC;                    % Total soil Capacity [m3]

SS_ratio=SSi/SSC;
DS_ratio=DSi/DSC;

if Tai > 0
    
    E = (0.051*(1-ALB_L)*Rsi*(Tai+9.5).^(1/2) - 2.4*(Rsi/Rai).^2 ...
             + 0.052*((Tai+20)*(1-RHi/100)*(1-0.38+0.54*WSi)))/1000;%m/day
    PET = (0.051*(1-ALB_G)*Rsi*(Tai+9.5).^(1/2)-2.4*(Rsi/Rai).^2 ...
        +0.048*(Tai+20)*(1-RHi/100)*(0.5+0.536*WSi))/1000; %m/day
    PET_demand = PET * CAi; petdemv(i)=PET_demand;
else
    E=0;
    PET=0; 
    PET_demand=0;
    petdemv(i)=PET_demand;
    
end

Fp = Pi/1000 * SAi * dt; fpv(i)=Fp;
Fsf = F_SF(Tai,Pi,CAi) * dt; fsfv(i)=Fsf;
Fr = F_R(Tai,Pi,CAi) * dt; frv(i)=Fr;
Fe= F_E(Tai,SAi,E) * dt;fev(i)=Fe;
Fsm = F_SM_daily(Tai,CAi,SPi,dt); fsmv(i)=Fsm;

Fslm = F_SLM_daily (SLi,monthi) * dt /30; fslmv(i)=Fslm;
Fdlm = F_DLM_daily (SVCi,SLi,monthi) * dt /30; fdlmv(i)=Fdlm;

Fdos = F_DOS_exp(DLi,HypsOutseep) * dt; fdosv(i)=Fdos;
Ftot = F_DOS_exp(LVi,HypsOutseep) * dt; 
Fsos = Ftot - Fdos; fsosv(i)=Fsos;

Fssd = F_SSD(SSi,SSC) * dt; fssdv(i)=Fssd;
Fdsd = F_DSD(DSi,DSC) * dt; fdsdv(i)=Fdsd;

Fsse = F_SSE_daily(SSi,PET_demand,dt) ; fssev(i)=Fsse;
Fdse = F_DSE_daily(DSi,PET_demand,Fsse,dt,TSC); fdsev(i)=Fdse;

%Fssi = F_SSI(Fr,Fsm,SSi,DSi,SSC,DSC) ;fssiv(i)=Fssi;
%Fro = F_RO(Fr,Fsm,SSi,DSi,SSC,DSC) ;frov(i)=Fro;

Frss = F_RSS(Fr,SS_ratio,DS_ratio); frssv(i)=Frss;
Frin = F_RIN(Fr,SS_ratio,DS_ratio); frinv(i)=Frin;
Fsmss = F_SMSS(Fsm,SS_ratio,DS_ratio); fsmssv(i)=Fsmss;
Fsmin = F_SMIN(Fsm,SS_ratio,DS_ratio); fsminv(i)=Fsmin;

Fin = F_IN(INi,C_IN) * dt; finv(i)=Fin;

Fof = F_OF(SLi,DLi,Max_V) * dt * Outflow_delay; fofv(i)=Fof;

if SSi + Frss + Fsmss < Fssd + Fsse
    Fssd = SSi * (Fssd/(Fssd+Fsse));
    Fsse = SSi * (Fsse/(Fssd+Fsse));
end

if DSi + Fssd < Fdsd + Fdse
    Fdsd = DSi * (Fdsd/(Fdsd+Fdse));
    Fdse = DSi * (Fdse/(Fdsd+Fdse));
end

if INi + Frin + Fsmin + Fdsd < Fin
    Fin = INi;
end

if SPi + Fsf < Fsmss + Fsmin
    Fdsd = SPi * (Fsmss/(Fsmss + Fsmin));
    Fdse = SPi * (Fsmin/(Fsmss + Fsmin));
end

%% Isotope equations

epsilon_sa=6.108*exp(17.27*Tai/(Tai+237.7));
epsilon_sw=6.108*exp(17.27*Twi/(Twi+237.7));
if RHi*(epsilon_sa/epsilon_sw)>=99
    h_n=99;
else 
    h_n=RHi*(epsilon_sa/epsilon_sw);
end

inverse_alpha=1/exp(0.35041*10^6/(Twi+273)^3-1.6664*10^3/(Twi+273)^2 ...
    +6.7123*(Twi+273)^-1-7.685*10^-3);

epsilon_eq=1000*(1-inverse_alpha);
epsilon_k=C*(1-h_n/100);
epsilon_tot=epsilon_eq+epsilon_k;

d_A=P_d18Oi-epsilon_eq;

limit_iso_enrichment=(RHi/100*d_A+epsilon_tot)/(RHi/100-10^(-3)*epsilon_tot);

d_L=SL_isoi/SLi;
dl(i)=d_L;
alpha_aragonite=exp((17.88*10^3/(Twi+273)-31.14)/10^3);
d18O_VSMOW_aragonite=alpha_aragonite*(1000+d_L)-1000;
d18O_VPDB_aragonite=(d18O_VSMOW_aragonite-30.92)/1.03092;
d18O_aragonite(i)=d18O_VPDB_aragonite;

if dayinyear>=summer_begin && dayinyear <summer_end
    i_summer=i_summer+1;
    
    temp_summer_arag(i_summer)=d18O_VPDB_aragonite;
    temp_summer_dl(i_summer)=d_L;
end

if dayinyear==summer_end && i_summer~=0
    year_summer_dl(i)=mean(temp_summer_dl);
    year_summer_arag(i)=mean(temp_summer_arag);
    temp_summer_arag=NaN(summer_len/dt,1);
    temp_summer_dl=NaN(summer_len/dt,1);
    i_summer=0;
end

if t==tf-1 && dayinyear <= summer_end
    year_summer_arag(i)=mean(temp_summer_arag(~isnan(temp_summer_arag)));
    year_summer_dl(i)=mean(temp_summer_dl(~isnan(temp_summer_dl)));
end

if d_L>=limit_iso_enrichment
    d_E=d_L;
elseif (inverse_alpha*d_L-h_n/100*d_A-epsilon_tot)/(1-h_n/100+0.001*epsilon_k)>=d_L
    d_E=d_L;
else
    d_E=(inverse_alpha*d_L-h_n/100*d_A-epsilon_tot)/(1-h_n/100+0.001*epsilon_k);
end

de(i)=d_E;

if DLi<1
    DL_d18O=0;
else 
    DL_d18O=DL_isoi/DLi;
end

if SSi<1
    SS_d18O=0;
else 
    SS_d18O=SS_isoi/SSi;
end

if DSi<1
    DS_d18O=0;
else 
    DS_d18O=DS_isoi/DSi;
end

if SPi<1
    SP_d18O=0;
else 
    SP_d18O=SP_isoi/SPi;
end

if INi<1
    IN_d18O=0;
else 
    IN_d18O=IN_isoi/INi;
end

%% Hydrologic mass balance equations for each reservoir

if MIX_OFF==true
    
    SLi = SLi + Fp + Fin - Fe  - Fsos - Fof - Fdos;
    SL(i+1)=SLi;
    
    SSi = SSi + Frss + Fsmss - Fsse - Fssd;
    SS(i+1)=SSi;
    
    DSi = DSi + Fssd - Fdse - Fdsd;
    DS(i+1)=DSi;

    INi = INi + Frin + Fsmin + Fdsd - Fin;
    IN(i+1)=INi;
    
    SPi = SPi + Fsf - Fsm;
    SP(i+1) = SPi;
    
    SL_isoi = SL_isoi + Fp*P_d18Oi + Fin*IN_d18O +  ...
        - Fe*d_E - d_L*(Fsos + Fof + Fdos);

    SL_iso(i+1)=SL_isoi;
    
    SS_isoi = SS_isoi + Frss*P_d18Oi + Fsmss*SP_d18O - Fsse*SS_d18O - Fssd*SS_d18O;
    SS_iso(i+1)=SS_isoi;

    DS_isoi = DS_isoi + Fssd*SS_d18O - Fdse*DS_d18O - Fdsd*DS_d18O;
    DS_iso(i+1)=DS_isoi;

    IN_isoi = IN_isoi + Frin*P_d18Oi + Fsmin*SP_d18O + Fdsd*DS_d18O - Fin*IN_d18O;
    IN_iso(i+1)=IN_isoi;

    SP_isoi = SP_isoi + Fsf*P_d18Oi - Fsm*SP_d18O;
    SP_iso(i+1) = SP_isoi;

else
    
    SLi = SLi + Fp + Fin + Fdlm - Fe - Fslm - Fsos - Fof;
    SL(i+1)=SLi;
    
    DLi = DLi + Fslm - Fdlm - Fdos;
    DL(i+1)=DLi;
    
    SSi = SSi + Frss + Fsmss - Fsse - Fssd;
    SS(i+1)=SSi;
    
    DSi = DSi + Fssd - Fdse - Fdsd;
    DS(i+1)=DSi;
 
    INi = INi + Frin + Fsmin + Fdsd - Fin;
    IN(i+1)=INi;
    
    SPi = SPi + Fsf - Fsm;
    SP(i+1) = SPi;

    % Isotope mass balance equations for each reservoir    

    SL_isoi = SL_isoi + Fp*P_d18Oi + Fin*IN_d18O + Fdlm*DL_d18O ...
        - Fe*d_E - d_L * (Fslm + Fsos + Fof);
    SL_iso(i+1)=SL_isoi;
    
    DL_isoi = DL_isoi + d_L*Fslm - Fdlm*DL_d18O - Fdos*DL_d18O;
    DL_iso(i+1)=DL_isoi;

    SS_isoi = SS_isoi + Frss*P_d18Oi + Fsmss*SP_d18O - Fsse*SS_d18O - Fssd*SS_d18O;
    SS_iso(i+1)=SS_isoi;

    DS_isoi = DS_isoi + Fssd*SS_d18O - Fdse*DS_d18O - Fdsd*DS_d18O;
    DS_iso(i+1)=DS_isoi;

    IN_isoi = IN_isoi + Frin*P_d18Oi + Fsmin*SP_d18O + Fdsd*DS_d18O - Fin*IN_d18O;
    IN_iso(i+1)=IN_isoi;

    SP_isoi = SP_isoi + Fsf*P_d18Oi - Fsm*SP_d18O;
    SP_iso(i+1) = SP_isoi;
end
%%
LVi=SLi+DLi;
Vv(i+1)=LVi;
i=i+1;

end
dates=datetime(years,months,days);
%% fix results and tidy up arrays

dates_clim=datetime(climate(:,3),climate(:,2),climate(:,1));
dates_clim_years=datetime(years(:),07,01);
dates_clim_years=unique(dates_clim_years);
dates_clim_years(isnat(dates_clim_years))=[];

daily_d18O_aragonite=d18O_aragonite(1:1/dt:end-1);

daily_dl=dl(1:1/dt:end-1);

daily_LDv=LDv(1:1/dt:end-1);
daily_LDv_mm=movmean(daily_LDv, [365/2 365/2]);

year_summer_arag(isnan(year_summer_arag))=[];
year_summer_arag_mm=movmean(year_summer_arag, [365/2 365/2]);

year_summer_dl(isnan(year_summer_dl))=[];
% winter starts end of sept
towinter=[9,30];
%summer starts end of march
tosummer=[3,31];

temp_winter_p=0.0;
temp_summer_p=0.0;
j=1;
%this function adds up precipitation over the whole dataset and divides it
%into summer and winter precipitation based on the dates above (yearly
%values)

for i=1:size(climate,1)
    if climate(i,2)<4 || climate(i,2)>=10
        temp_winter_p=temp_winter_p+climate(i,4);
    else
        temp_summer_p=temp_summer_p+climate(i,4);
    end
    if [climate(i,2),climate(i,1)] == towinter
        year_summer_p(j)=temp_summer_p;
        temp_summer_p=0;
        j=j+1;
    elseif [climate(i,2),climate(i,1)] == tosummer
        year_winter_p(j)=temp_winter_p;
        temp_winter_p=0;
    end
    if i==size(climate,1)
        year_summer_p(end)=temp_summer_p;
    end
end
%%
%load level logger lake level
LLobs=readmatrix('castorLakeLevel0519.csv');
dates_obs=datetime(LLobs(:,3),LLobs(:,2),LLobs(:,1));
LLobs_mm=movmean(LLobs(:,4), [365/2 365/2]);

%load satellite image lake level
LLsat=readmatrix('SatDepths.csv');
dates_sat=datetime(LLsat(:,3),LLsat(:,2),LLsat(:,1));

%load 18O sample values
d18O_samples=readmatrix('18O samples Castor 2.csv');
dates_samples=datetime(d18O_samples(:,3),d18O_samples(:,2),d18O_samples(:,1));

%load Sam's annual 100-year 18O record
d18O_core=readmatrix('18OannualSam.csv');
dates_core=datetime(d18O_core(:,1),07,01);

%% linear models for determining error

% calculate rmse for level logger observations
[~,a1,b1]=intersect(dates_clim,dates_obs);
lm_obs=fitlm(daily_LDv_mm(a1),LLobs_mm(b1));
obs_rmse=rmse(daily_LDv_mm(a1),LLobs_mm(b1));
obs_corr=corr(daily_LDv_mm(a1),LLobs_mm(b1));
obs_mean_diff=mean(daily_LDv_mm(a1))-mean(LLobs_mm(b1));

%calculate rmse for 18O water samples
[~,a1,b1]=intersect(dates_clim,dates_samples);
% lm_samples=fitlm(d18O_samples(b1,4),daily_dl(a1));
samples_rmse=rmse(d18O_samples(b1,4),daily_dl(a1));

%calculate rmse for 18O sediment values 
trail_avg_year_summer_arag=movmean(year_summer_arag,[10 0]);
[~,a1,b1]=intersect(dates_clim_years,dates_core);
core_rmse=rmse(trail_avg_year_summer_arag(a1),movmean(d18O_core(b1,2),[10 0]));
core_corr=corr(trail_avg_year_summer_arag(a1),movmean(d18O_core(b1,2),[10 0]));

if length(a1) > length(trail_avg_year_summer_arag)
    bad_run=true;
else
    lm_core=fitlm(d18O_core(b1,2),trail_avg_year_summer_arag(a1));
    bad_run=false;
end


meanCore_meanYearArag_diff=mean(year_summer_arag(a1))-mean(d18O_core(b1,2));

% plotResults20thCentSingleRun;
% %calculate amounts of water lost through evaporation, outseepage, overflow,
% %water lost total and proportions with respect to total for each
% F_evap=sum(fev(~isnan(fev)));
% F_os=sum(fsosv(~isnan(fsosv)))+sum(fdosv(~isnan(fdosv)));
% F_of=sum(fofv(~isnan(fofv)));
% F_out=F_evap+F_os+F_of;
% evap_ratio=F_evap/F_out;
% os_ratio=F_os/F_out;
% of_ratio=F_of/F_out;
