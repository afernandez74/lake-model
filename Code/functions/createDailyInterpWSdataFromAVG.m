function dailyAVGWS = createDailyInterpWSdataFromAVG()

load ('../Data/AVGclimate.mat','WS')

Y=NaN(366,1);
% Y(:,1)=1:365;
Y(1)=WS(1);%jan
Y(32)=WS(2);%feb
Y(60)=WS(3);%mar
Y(91)=WS(4);%abr
Y(121)=WS(5);%may
Y(152)=WS(6);%jun
Y(182)=WS(7);%jul
Y(213)=WS(8);%agu
Y(244)=WS(9);%sep
Y(274)=WS(10);%oct
Y(305)=WS(11);%nov
Y(335)=WS(12);%dec
Y(size(Y,1))=Y(1);
t=[1:366]';

Y2=Y;
t2=t;

Y2(isnan(Y)==1)=[];
t2(isnan(Y)==1)=[];

dailyAVGWS=interp1(t2,Y2,t);
dailyAVGWS(size(dailyAVGWS,1))=[];
temp=dailyAVGWS(size(dailyAVGWS,1)-13:size(dailyAVGWS,1));
dailyAVGWS(size(dailyAVGWS,1)-13:size(dailyAVGWS,1))=[];
dailyAVGWS=[temp;dailyAVGWS];
end
