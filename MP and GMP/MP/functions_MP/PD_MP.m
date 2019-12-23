function x_opt = PD_MP(y, coef_Matrix, memory_deg, memory_depth)

coef_vec  = reshape(coef_Matrix, memory_deg*(memory_depth+1), 1);
input_len = length(y);
first_n   = memory_depth+1;
x_opt     = zeros(input_len-(memory_depth+1),1);

for n = first_n:input_len
    current_time = zeros(1, memory_deg*(memory_depth+1));
    j            = 1;
    for m = 0:memory_depth
        for k = 1:memory_deg
            current_time(j) =  y(n-m)*(abs(y(n-m))^(k-1));
            j               = j+1;
        end
    end
    x_opt(n-first_n+1) = current_time*coef_vec;
end

end