%creates yearly mixing depth data with a max mixing depth 

Y=NaN(366,1);
Y(1)=10;%jan
Y(32)=10;%feb
Y(60)=10;%mar
if MD_max*0.083 <= 10
    Y(91)=10;
else
    Y(91)=MD_max*0.083;%abr
end
if MD_max*0.25 <= 10
    Y(121)=10;
else
    Y(121)=MD_max*0.25;%may
end
if MD_max*0.5 <= 10
    Y(152)=10;
else
    Y(152)=MD_max*0.5;%jun
end

Y(182)=MD_max;%jul
Y(213)=MD_max;%agu
Y(244)=MD_max;%sep
Y(274)=MD_max;%oct
Y(305)=10;%nov
Y(335)=10;%dec
Y(size(Y,1))=Y(1);
t=[1:366]';

Y2=Y;
t2=t;

Y2(isnan(Y)==1)=[];
t2(isnan(Y)==1)=[];

dailyAVGMD=interp1(t2,Y2,t);
dailyAVGMD(size(dailyAVGMD,1))=[];
temp=dailyAVGMD(size(dailyAVGMD,1)-13:size(dailyAVGMD,1));
dailyAVGMD(size(dailyAVGMD,1)-13:size(dailyAVGMD,1))=[];
dailyAVGMD=[temp;dailyAVGMD];

