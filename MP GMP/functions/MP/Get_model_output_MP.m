function Amp_out = Get_model_output_MP(coef_Matrix, input, memory_deg, memory_depth)

coef_vec   = reshape(coef_Matrix, memory_deg*(memory_depth+1), 1);
input_len  = length(input);
first_n    = memory_depth+1;
Amp_out    = zeros(input_len-memory_depth,1);

for n = first_n:input_len
    current_time = zeros(1, memory_deg*(memory_depth+1));
    j            = 1;
    for m = 0:memory_depth
        for k = 1:memory_deg
            current_time(j) =  input(n-m)*(abs(input(n-m))^(k-1));
            j               = j+1;
        end
    end
    Amp_out(n-first_n+1) = current_time*coef_vec;
end

end