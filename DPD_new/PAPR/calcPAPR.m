function x_papr = calcPAPR(x)

    avg = 0;
    peak = 0;
    for i=1:length(x)
        curr = abs(x(i));
        if peak<curr
            peak = curr;
        end
    end
    papr = (peak^2)/(rms(x)^2);
    papr_db = 10*log10(papr);
    x_papr = papr_db;
end
