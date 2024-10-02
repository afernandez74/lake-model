function F_DSD = F_DSD (RES_DS,DSC)

% Deep soil drainage m3

    if RES_DS - DSC < 0

        F_DSD = 0;
    else

        F_DSD = RES_DS - DSC ;

    end

end

