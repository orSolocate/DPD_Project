function output_err_vec = Get_output_err_vec(y_d, y)

    delay           = finddelay(y, y_d);
    if(delay >=0)
        output_err_vec = [zeros(delay,1); y_d(delay+1:end) - y(1:end-delay)];
    end
    if(delay < 0)
        output_err_vec = [zeros(-delay,1); y_d(1:end+delay) - y(1-delay:end)];
    end
    
end