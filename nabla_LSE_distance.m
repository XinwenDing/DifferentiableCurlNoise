function del_lse = nabla_LSE_distance(vx, vy, x, y)
% This function calcuates the gradient of the LSE distance function.
% vx and vy are the x- and y-coordinates of polygon vertices.
% x and y are location of query point.
    [dx, dy] = deal(1e-4, 1e-4);
    [d_pdx, ~] = min_poly_dist(vx, vy, x + dx, y);
    [d_mdx, ~] = min_poly_dist(vx, vy, x - dx, y);
    [d_pdy, ~] = min_poly_dist(vx, vy, x, y + dy);
    [d_mdy, ~] = min_poly_dist(vx, vy, x, y - dy);
    del_x = (d_pdx - d_mdx) / (2 * dx);
    del_y = (d_pdy - d_mdy) / (2 * dy);
    del_lse = [del_x, del_y];
end
