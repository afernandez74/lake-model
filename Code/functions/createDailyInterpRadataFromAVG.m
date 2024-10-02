function dailyAVGRa = createDailyInterpRadataFromAVG()

load ('../Data/AVGclimate.mat','Ra')

Y=NaN(366,1);
% Y(:,1)=1:365;
Y(1)=Ra(1);%jan
Y(32)=Ra(2);%feb
Y(60)=Ra(3);%mar
Y(91)=Ra(4);%abr
Y(121)=Ra(5);%may
Y(152)=Ra(6);%jun
Y(182)=Ra(7);%jul
Y(213)=Ra(8);%agu
Y(244)=Ra(9);%sep
Y(274)=Ra(10);%oct
Y(305)=Ra(11);%nov
Y(335)=Ra(12);%dec
Y(size(Y,1))=Y(1);
t=[1:366]';

Y2=Y;
t2=t;

Y2(isnan(Y)==1)=[];
t2(isnan(Y)==1)=[];

dailyAVGRa=interp1(t2,Y2,t);
dailyAVGRa(size(dailyAVGRa,1))=[];
temp=dailyAVGRa(size(dailyAVGRa,1)-13:size(dailyAVGRa,1));
dailyAVGRa(size(dailyAVGRa,1)-13:size(dailyAVGRa,1))=[];
dailyAVGRa=[temp;dailyAVGRa];
end
