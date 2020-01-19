function err_vec = Get_err_vec(x, y)

    delay             = finddelay(x, y);
    if(delay >=0)
        err_vec  = [zeros(delay,1);((y(delay+1:end)) - x(1:end-delay))];
    end
    if(delay < 0)
        err_vec  = [((y(1:end+delay)) - x(1-delay:end));zeros(delay,1)];
    end

end