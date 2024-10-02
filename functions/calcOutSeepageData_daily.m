function HypsOS = calcOutSeepageData_daily (a,b,Hyps)
%Script that calculates outseepage rates for basin

%This script takes the hypsometry file (with depth, surface area and
%volume), and appends columns of data shown below:
% 1          2                  3           4              5        
% Depth [cm] Surface Area [m2]	Volume [m3] Est Radius [m] Dist from shore [m]
% 6                 7                         8                    9
% npSA(diff)[m2]    Outseepage Rate[cm/day]   Outseepage [m3/day]  Accumulated Outseepage w/depth [m3]

%with these calculations, the outseepage value for a determined lake depth
%or volume can be interpolated from the table therefore just calculating
%all outseepage values just once, and using a table look up function for
%the model.

[r,c]=size(Hyps);
A=Hyps;
A(r,c+5)=0;
A(:,4)=(A(:,2)/pi()).^(0.5); %radius [m]
A(:,5)=flipud(A(:,4)); %distance from shore [m]

dH=(A(2,1)-A(3,1))/100;
npSA=NaN(r,1);
npSA(r)=0;
Outseepage=NaN(r,1);
Outseepage(r)=0;

% Outseepage rate at each radius
% OR=a*exp(Dist from shore/b) 

% HypsOutseep(:,5)=A*exp(HypsOutseep(:,4)/b);

for i=r-1:-1:1
    if i==1
        dH=(A(1,1)-A(2,1))/100;
    end
    r_1=A(i,4); % radius of i-th slice
    r_2=A(i+1,4); % radius of previous slice
    r_diff=r_1-r_2; % radius difference
    l=(r_diff^2+dH^2)^(0.5); % non planar width
    npSA(i)=pi()*(r_1+r_2)*l; % array of non planar areas for each radius
    if a*exp(-(A(i,5))/b)>100.0 % exponential equation for outseepage flux
       Outseepage(i)= 100.0;
    elseif a*exp(-(A(i,5))/b)<0.01
       Outseepage(i)= 0.01;
    else
       Outseepage(i)=a*exp(-(A(i,5))/b);
    end
    A(i,6)=npSA(i);
    
end

A(:,7)=Outseepage;
A(:,8)=(A(:,7)/100).*(A(:,6));
A(r,9)=A(r,8);
A(r,10)=r;

for i=r-1:-1:1
    A(i,9)=A(i+1,9)+A(i,8);
    A(i,10)=i;
end

HypsOS=A;
end