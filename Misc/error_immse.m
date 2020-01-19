function error_vector = error_immse(y_d,y_no_PD,output_MP,output_GMP, memory_depth)

error_vector=zeros(3,length(y_no_PD),'double')
y_d_for_error=y_d(memory_depth+1:end)
output_MP_for_error=output_MP(memory_depth+1:end)
output_GMP_for_error=output_GMP(memory_depth+1:end)

for i=1:3
    for j=1:length(y_no_PD)
        switch(i)
            case 1
                error_vector(i,j)=immse(y_no_PD(1:j),y_d_for_error(1:j));
                
            case 2
                error_vector(i,j)=immse(output_MP_for_error(1:j),y_d_for_error(1:j));
            case 3
                error_vector(i,j)=immse(output_GMP_for_error(1:j),y_d_for_error(1:j));
                return
        end
    end
end
end

