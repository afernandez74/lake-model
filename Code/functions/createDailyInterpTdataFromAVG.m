function dailyAVGTa = createDailyInterpTdataFromAVG()

load ('../Data/AVGclimate.mat','Ta')

Y=NaN(366,1);
% Y(:,1)=1:365;
Y(1)=Ta(1);%jan
Y(32)=Ta(2);%feb
Y(60)=Ta(3);%mar
Y(91)=Ta(4);%abr
Y(121)=Ta(5);%may
Y(152)=Ta(6);%jun
Y(182)=Ta(7);%jul
Y(213)=Ta(8);%agu
Y(244)=Ta(9);%sep
Y(274)=Ta(10);%oct
Y(305)=Ta(11);%nov
Y(335)=Ta(12);%dec
Y(size(Y,1))=Y(1);
t=[1:366]';

Y2=Y;
t2=t;

Y2(isnan(Y)==1)=[];
t2(isnan(Y)==1)=[];

dailyAVGTa=interp1(t2,Y2,t);
dailyAVGTa(size(dailyAVGTa,1))=[];
temp=dailyAVGTa(size(dailyAVGTa,1)-13:size(dailyAVGTa,1));
dailyAVGTa(size(dailyAVGTa,1)-13:size(dailyAVGTa,1))=[];
dailyAVGTa=[temp;dailyAVGTa];
end
