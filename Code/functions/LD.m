function LD = LD (Lake_Vol_m3,hyp)

% Lake Depth cm
% Calculates lake depth from Hypsometry data with lake volume [cm]

% If lake volume is greater than the greatest value in the hypsography
% spreadsheet, returns the greatest depth value 

    if Lake_Vol_m3 <= hyp(1,3)
        
        LD = interp1(hyp(:,3),hyp(:,1),Lake_Vol_m3);

    else
        
        LD = hyp(1,1);
        
    end
        
end
