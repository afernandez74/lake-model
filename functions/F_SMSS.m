function F_SMSS = F_SMSS(Fsm,SS_ratio,DS_ratio)

if SS_ratio < 1 
    
    F_SMSS=Fsm;
    
elseif SS_ratio >= 1 && DS_ratio < 1
    
    F_SMSS=Fsm/2;

else
    
    F_SMSS=0;
end
end
