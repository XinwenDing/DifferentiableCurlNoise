function best_pos = adhere_isocontour(init_p, psi, d0, curr_pos, obj_list, center, vertex, cr, f, g)
% With our construction, the potential psi' is differentiable, which allows
% us to force particles to locate on the same isocontour all the time using
% a first-order variation of Newton's Method (Chopp [2001]).
% init_p: initial potential value
% psi: psi' in our paper
% d0: self-adaptive influence radius
% curr_pos: current query point position. size(curr_pos) = particle_num x 2. 
%           The first column of curr_pos stores the x-coordinates of the query 
%           points, and the second column of curr_pos stores the y-coordinates of 
%           the query points.
% obj_list: list of objects.
% center: e.g. center = {cx, cy}; center of the objects, where cx and cy
%         stores the x- and y-coordinates of the center of the objects
% vertex: e.g. vertex = {vx, vy}; vertices of the objects, where vx and vy
%         stores the x- and y-coordinates of the vertices of the objects
% cr: Radius of objects, if properly defined (e.g. circles and regular polygons have radius)
% f: SC-exterior map
% g: inverse of f.
% For definition and data structure of obj_list, center, vertex, cr, please
% search (Ctrl-F) and check the "Obstacle info" part in any file starts with "flow particle"
% Please check out our paper for more details.
% WARNING: This function is EXTREMELY expensive.
    gradient = nabla_particle_potential(psi, d0, curr_pos, obj_list, center, vertex, cr, f, g);
    norm_square = vecnorm(gradient, 2, 2).^2;
    p = get_particle_potential(psi, d0, curr_pos, obj_list, center, vertex, cr, f, g);
    diff_potential = p - init_p;
    k = 1;
    diff_min = 100;
    best_pos = curr_pos;
    while max(abs(diff_potential)) > 1e-5 && k <= 10
        curr_pos = curr_pos - diff_potential .* gradient ./ norm_square;
        gradient = nabla_particle_potential(psi, d0, curr_pos, obj_list, center, vertex, cr, f, g);
        norm_square = vecnorm(gradient, 2, 2).^2;
        p = get_particle_potential(psi, d0, curr_pos, obj_list, center, vertex, cr, f, g);
        diff_potential = p - init_p;
        if max(abs(diff_potential)) < diff_min
            best_pos = curr_pos;
            diff_min = max(abs(diff_potential));
        end
        k = k + 1;
    end
end