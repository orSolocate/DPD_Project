function Amp_out = Get_model_output_GMP(coef_Vector, input, orders)

input_len       = length(input);
first_n         = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);
last_n          = input_len-orders(8);
columns_aligned = (1+orders(1))*orders(2);
columns_lagging = (1+orders(3))*(orders(4)-1)*orders(5);
columns_leading = (1+orders(6))*(orders(7)-1)*orders(8);
Amp_out         = zeros(last_n-first_n,1);

for n = first_n:last_n
    current_time = zeros(1, columns_aligned+columns_lagging+columns_leading);
    j            = 1;
    for m = 0:orders(1)
        for k = 1:orders(2)
            current_time(n-first_n+1,j) =  input(n-m)*(abs(input(n-m))^(k-1));
            j                          = j+1;
        end
    end
 
    for m = 0:orders(3)
        for k = 2:orders(4)
            for p=1:orders(5)
                current_time(n-first_n+1,j) =  input(n-m)*(abs(input(n-m-p))^(k-1));
                j                          = j+1;
            end
        end
    end
   
    for m = 0:orders(6)
        for k = 2:orders(7)
            for q=1:orders(8)
                current_time(n-first_n+1,j) =  input(n-m)*(abs(input(n-m+q))^(k-1));
                j                          = j+1;
            end
        end
    end  
    
    Amp_out(n-first_n+1) = current_time*coef_Vector;
end

end