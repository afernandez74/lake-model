% calculate the corresponding day number within the year(1-365) 
% for a month and day input

function dayofyear = day365(month,day)
temp=0;
if month==1
    temp=0;
elseif month <= 7 
    for i = 1 : month-1
        if mod(i,2)==1
            temp=temp+31;
        elseif mod(i,2)==0 && i == 2
            temp=temp+28;
        else
            temp=temp+30;
        end
    end
else
    temp=212;
    for i = 8 : month-1
        if mod(i,2)==1
            temp=temp+30;
        else
            temp=temp+31;
        end
    end
end
    dayofyear=temp+day;
end

