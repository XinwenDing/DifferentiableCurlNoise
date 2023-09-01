function [ramp_input, obj_idx, cp] = potential_prep(d0, ptx, pty, obj_list, center, vertex, cr)
% This function do the prep work for the function "get_particle_potential".
% The function takes a series of query point, where x-coordinates are stored in ptx, and y-coordinates are
% stored in pty, and the constant influence radius d0
% For each query point, the function outputs the value of the ramp function, the closest oject, 
% and the closest point (cp) on the boundry of the object
% d0: constant influence radius
% ptx: x-coordinates of a set of query points
% pty: y-coordinates of a set of query points
% obj_list: list of objects.
% center: e.g. center = {cx, cy}; center of the objects, where cx and cy
%         stores the x- and y-coordinates of the center of the objects
% vertex: e.g. vertex = {vx, vy}; vertices of the objects, where vx and vy
%         stores the x- and y-coordinates of the vertices of the objects
% cr: Radius of objects, if properly defined (e.g. circles and regular polygons have radius)
% For definition and data structure of obj_list, center, vertex, cr, please
% search (Ctrl-F) and check the "Obstacle info" part in any file starts with "flow particle"
    obj_num = length(obj_list);
    cx = center{1}; % eg. cx = [x1, x2];
    cy = center{2}; % eg. cy = [y1, y2];
    vx = vertex{1}; % eg. vx = {"None", [vx2]};
    vy = vertex{2}; % eg. vy = {"None", [vy2]};
    closest_dist = zeros(obj_num, 1);
    loc = zeros(obj_num, 2);
    for i = 1 : obj_num
        if obj_list(i) == "circle"
            [closest_dist(i), loc(i,:)] = min_circ_dist(cx(i), cy(i), cr(i), ptx, pty);
        elseif obj_list(i) == "polygon"
            [closest_dist(i), loc(i,:), ~] = min_poly_dist(vx{i}, vy{i}, ptx, pty);
        end
    end
    [~, obj_idx] = min(closest_dist);
    if obj_list(obj_idx) == "circle" || closest_dist(obj_idx) < 0
        cp = loc(obj_idx,:);
    else
        cp = NaN;
    end
    
    a1 = 200; 
    a2 = 200;
    smooth_min_below = @(x) log(sum(exp(-a1 .* x))) / (-a1);
    smooth_min_above = @(x) sum(x .* exp(-a1 * x)) / sum(exp(-a1 * x));
    smooth_min2sum = @(x) (length(x) - 2) * log(sum(exp(-a2 .* x))) / a2 - ...
                            sum(log(sum(exp(-a2 * x)) - exp(-a2 * x))) / a2;
    approx_min = smooth_min_above(closest_dist);
    approx_min2sum = smooth_min2sum(closest_dist);
    approx_d2 = approx_min2sum - approx_min;
    act_d0 = smooth_min_below([approx_d2, d0]);
    ramp_input = approx_min / act_d0;
end