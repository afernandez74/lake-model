function F_RSS = F_RSS(Fr,SS_ratio,DS_ratio)

if SS_ratio < 1
    
    F_RSS=Fr;
    
elseif SS_ratio >= 1 && DS_ratio <1
    
    F_RSS=Fr/2;

else
    
    F_RSS=0;
end

end
