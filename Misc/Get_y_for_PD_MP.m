function y_for_PD  = Get_y_for_PD_MP(y_d, output, position, memory_depth)

y_for_PD(memory_depth+1) = y_d(position);
for n=1:memory_depth
    y_for_PD(n) = output(position-memory_depth-1+n);
end

end