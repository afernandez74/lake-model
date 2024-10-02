function F_DOS_exp = F_DOS_exp (DLV,hyp)

% Deep lake outseepage m3/month
% Using exponential outseepage function

if DLV >= hyp(end,3)

    if DLV <= hyp(1,3) 
            F_DOS_exp = interp1(hyp(:,3),hyp(:,9),DLV);
    else 
            F_DOS_exp = hyp(1,9);
    end
    
    else
        
        F_DOS_exp =0;
    
end
end

