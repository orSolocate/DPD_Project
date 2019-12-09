function updated_coef = Update_coef_GMP(x, y, PD_coef_Vec_GMP_WL, orders_WL, miu)

    delay           = finddelay(x, y);
    if(delay >=0)
        updated_coef = PD_coef_Vec_GMP_WL +...
            miu.*Get_coef_GMP(y(delay+1:end)', x(1:end-delay)', orders_WL);
    end
    if(delay < 0)
        updated_coef = PD_coef_Vec_GMP_WL +...
            miu.*Get_coef_GMP(y(1:end+delay)', x(1-delay:end)', orders_WL);
    end
    
end