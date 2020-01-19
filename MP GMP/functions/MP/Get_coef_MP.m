function coef_Matrix = Get_coef_MP(input, output, memory_deg, memory_depth)

input_len  = length(input);
first_n    = memory_depth+1;
phi_matrix = zeros(input_len-first_n, memory_deg*(memory_depth+1)); 

for n = first_n:input_len
    j = 1; 
    for m = 0:memory_depth
        for k = 1:memory_deg
            phi_matrix(n-first_n+1,j) =  input(n-m)*(abs(input(n-m))^(k-1));
            j                         = j+1;
        end
    end
end

coef_vec    = phi_matrix\(output(first_n:input_len).');
coef_Matrix = reshape(coef_vec,memory_depth+1,memory_deg);

end
