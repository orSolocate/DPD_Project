%orders is a vector in the length of 8 with all orders in the sums of the GMP
%orders(1) - Ma: memory depth of the aligned polynomial
%orders(2) - Ka: nonlinearity order of the aligned polynomial
%orders(3) - Mb: memory depth of the lagging polynomial
%orders(4) - Kb: nonlinearity order of the lagging polynomial
%orders(5) - P : cross-terms order of the lagging polynomial
%orders(6) - Mc: memory depth of the leading polynomial
%orders(7) - Kc: nonlinearity order of the leading polynomial
%orders(8) - Q : cross-terms order of the leading polynomial

function coef_Vector = Get_coef_GMP(input, output, orders)

columns_aligned = (1+orders(1))*orders(2);
columns_lagging = (1+orders(3))*(orders(4)-1)*orders(5);
columns_leading = (1+orders(6))*(orders(7)-1)*orders(8);
input_len       = length(input);
first_n         = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);
last_n          = input_len-orders(8);
rows            = 1+last_n-first_n;
phi_aligned     = zeros(rows, columns_aligned);
phi_lagging     = zeros(rows, columns_lagging);
phi_leading     = zeros(rows, columns_leading);

for n = first_n:last_n
    j = 1; 
    for m = 0:orders(1)
        for k = 1:orders(2)
            phi_aligned(n-first_n+1,j) =  input(n-m)*(abs(input(n-m))^(k-1));
            j                          = j+1;
        end
    end
end

for n = first_n:last_n
    j = 1; 
    for m = 0:orders(3)
        for k = 2:orders(4)
            for p=1:orders(5)
                phi_lagging(n-first_n+1,j) =  input(n-m)*(abs(input(n-m-p))^(k-1));
                j                          = j+1;
            end
        end
    end
end

for n = first_n:last_n
    j = 1; 
    for m = 0:orders(6)
        for k = 2:orders(7)
            for q=1:orders(8)
                phi_leading(n-first_n+1,j) =  input(n-m)*(abs(input(n-m+q))^(k-1));
                j                          = j+1;
            end
        end
    end
end

phi_all = [phi_aligned, phi_lagging, phi_leading];
coef_Vector    = phi_all\(output(first_n:last_n).');

end
