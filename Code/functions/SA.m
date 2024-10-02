function SA = SA (Lake_Vol_m3,hyp)

% Lake Surface Area m2
% Calculates Surface area from Hypsometry data with lake volume input [m2]

% If lake volume is greater than the greatest value in the hypsography
% spreadsheet, returns the greatest surface area value 

    if Lake_Vol_m3 <= hyp(1,3)
        
        SA = interp1(hyp(:,3),hyp(:,2),Lake_Vol_m3);

    else
        
        SA = hyp(1,2);
        
    end
        
end
