function updated_coef = Update_coef_MP(y, x, PD_coef_Matrix_WL, mem_deg, mem_depth, miu)

    delay             = finddelay(x, y);
    if(delay >=0)
        u = x(1:end-delay);
        v = y(delay+1:end);
    end
    if(delay < 0)
        u = x(1-delay:end);
        v = y(1:end+delay);
    end

    input_len  = length(v);
    first_n    = mem_depth+1;
    Y          = zeros(input_len-first_n, mem_deg*(mem_depth+1)); 

    for n = first_n:input_len
        j = 1; 
        for m = 0:mem_depth
            for k = 1:mem_deg
                Y(n-first_n+1,j) =  v(n-m)*(abs(v(n-m))^(k-1));
                j                = j+1;
            end
        end
    end
    
    coef_vec     = Y\(u(first_n:input_len));
    %coef_vec     = inv(ctranspose(Y)*Y)*(ctranspose(Y)*u(first_n:input_len));
    coef_Matrix  = reshape(coef_vec,mem_depth+1,mem_deg);
    updated_coef = PD_coef_Matrix_WL + miu.*coef_Matrix;
    
%    delay             = finddelay(x, y);
%    if(delay >=0)
%        updated_coef  = PD_coef_Matrix_WL +...
%            miu.*Get_coef_MP(x(delay+1:end)', y(1:end-delay)', mem_deg, mem_depth);
%    end
%    if(delay < 0)
%        updated_coef  = PD_coef_Matrix_WL +...
%            miu.*Get_coef_MP(x(1:end+delay)', y(1-delay:end)', mem_deg, mem_depth);
%    end

end