function dailyAVGMD = createDailyInterpMDdataFromAVG()

load ('../Data/AVGclimate.mat','MD')

Y=NaN(366,1);
% Y(:,1)=1:365;
Y(1)=MD(1);%jan
Y(32)=MD(2);%feb
Y(60)=MD(3);%mar
Y(91)=MD(4);%abr
Y(121)=MD(5);%may
Y(152)=MD(6);%jun
Y(182)=MD(7);%jul
Y(213)=MD(8);%agu
Y(244)=MD(9);%sep
Y(274)=MD(10);%oct
Y(305)=MD(11);%nov
Y(335)=MD(12);%dec
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
end
