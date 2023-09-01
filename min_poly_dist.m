function [min_dist, cp, flip] = min_poly_dist(vx, vy, x, y)
% This function finds the min distance and the closest point (cp) from the 
% query point (x,y) to a polygon, where vx stores the x-coordinates of its 
% vertices and vy stores the y-coordinates of its vertices
% flip is a variable for the tweek particle function, indicating the
% direction of pushing the particle.
    [min_dist, pos] = LSE_distance(vx, vy, x, y);
    cp = pos;
    flip = 0;
    [in, on] = inpolygon(x, y, vx, vy);
    if in == 1 && on == 0
        min_dist = -abs(min_dist);
        flip = 1;
    end
end