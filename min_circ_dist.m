function [dist, cp] = min_circ_dist(x0, y0, R, x, y)
% This function finds the min distance, and the closest point (cp)
% from the query point (x,y) to a circle centered at (x0, y0) with radius R.
    vector_to_center = [(x - x0) (y - y0)];
    dist_to_circle = @(x, y) norm(vector_to_center) - R;
    dist = dist_to_circle(x,y);
    cp = [x0 y0] + R * (vector_to_center / norm(vector_to_center));
end