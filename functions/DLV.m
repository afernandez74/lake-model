function DLV = DLV (dld,hyp)

% Deep Lake Volume m3
% Calculates Deep Lake Volume from Hypsometry data with Deep Lake Depth [m3]

% If Deep Lake Depth is greater than the greatest value in the hypsography
% spreadsheet, returns the greatest depth value
% If Deep Lake Depth value is negative, returns 0

    if dld <= hyp(1,1) && dld > 0
        
        DLV = interp1(hyp(:,1),hyp(:,3),dld);

    elseif dld <= 0
        
        DLV = 0;
        
    else
            
        DLV = hyp(1,3);
        
    end
        
end
