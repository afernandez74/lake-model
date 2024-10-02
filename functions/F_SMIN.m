function F_SMIN = F_SMIN(Fsm,SS_ratio,DS_ratio)

if SS_ratio >= 1 && DS_ratio >= 1 
    
    F_SMIN=Fsm;
    
elseif SS_ratio >= 1 && DS_ratio < 1
    
    F_SMIN=Fsm/2;

else
    
    F_SMIN=0;
end
end
