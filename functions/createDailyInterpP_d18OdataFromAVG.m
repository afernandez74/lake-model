function dailyAVGP_d18O = createDailyInterpP_d18OdataFromAVG()

load ('../Data/AVGclimate.mat','P_d18O')

Y=NaN(366,1);
% Y(:,1)=1:365;
Y(1)=P_d18O(1);%jan
Y(32)=P_d18O(2);%feb
Y(60)=P_d18O(3);%mar
Y(91)=P_d18O(4);%abr
Y(121)=P_d18O(5);%may
Y(152)=P_d18O(6);%jun
Y(182)=P_d18O(7);%jul
Y(213)=P_d18O(8);%agu
Y(244)=P_d18O(9);%sep
Y(274)=P_d18O(10);%oct
Y(305)=P_d18O(11);%nov
Y(335)=P_d18O(12);%dec
Y(size(Y,1))=Y(1);
t=[1:366]';

Y2=Y;
t2=t;

Y2(isnan(Y)==1)=[];
t2(isnan(Y)==1)=[];

dailyAVGP_d18O=interp1(t2,Y2,t);
dailyAVGP_d18O(size(dailyAVGP_d18O,1))=[];
temp=dailyAVGP_d18O(size(dailyAVGP_d18O,1)-13:size(dailyAVGP_d18O,1));
dailyAVGP_d18O(size(dailyAVGP_d18O,1)-13:size(dailyAVGP_d18O,1))=[];
dailyAVGP_d18O=[temp;dailyAVGP_d18O];
end
