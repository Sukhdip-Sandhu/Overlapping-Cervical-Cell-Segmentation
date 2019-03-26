function [color] = edge_LUT(i)
num_colors = 3;
r = mod(i,num_colors);

if r == 0
    color = 'r';
    return
end

if r == 1
    color = 'g';
    return
end

if r == 2
    color = 'b';
    return
end

end

