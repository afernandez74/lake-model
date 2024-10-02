% calculate the corresponding day number within the year(1-365) 
% for a month and day input

function dayofyear = day365(month,day)

    m=num2str(month);
    d=num2str(day);
    date=strcat(m,'-',d,'-','1');
    dayofyear=days365('1-1-1',date)+1;
end
