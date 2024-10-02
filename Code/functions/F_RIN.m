function F_RIN = F_RIN(Fr,SS_ratio,DS_ratio)

if SS_ratio >= 1 && DS_ratio >= 1
    
    F_RIN=Fr;
    
elseif SS_ratio >= 1 && DS_ratio < 1
    
    F_RIN=Fr/2;

else
    
    F_RIN=0;
end
end
