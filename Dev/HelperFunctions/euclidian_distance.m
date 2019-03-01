%{
    Function that returns Euclidian Distance between two (X,Y) coordinates.
%}
function [dist] = euclidian_distance(p1, p2)
diff = p2 - p1;
dist = sqrt(diff * diff');
end

