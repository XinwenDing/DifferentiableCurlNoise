function gradient = nabla_particle_potential(psi, d0, curr_pos, obj_list, center, vertex, cr, f, g)
% This function calcuates the gradient of the psi (psi' in paper).
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
    dx = 1e-4;
    dy = 1e-4;
    pos = [curr_pos(:,1) + dx, curr_pos(:,2)];
    p_pdx = get_particle_potential(psi, d0, pos, obj_list, center, vertex, cr, f, g);
    pos = [curr_pos(:,1) - dx, curr_pos(:,2)];
    p_mdx = get_particle_potential(psi, d0, pos, obj_list, center, vertex, cr, f, g);
    pos = [curr_pos(:,1), curr_pos(:,2) + dy];
    p_pdy = get_particle_potential(psi, d0, pos, obj_list, center, vertex, cr, f, g);
    pos = [curr_pos(:,1), curr_pos(:,2) - dy];
    p_mdy = get_particle_potential(psi, d0, pos, obj_list, center, vertex, cr, f, g);
    del_x = (p_pdx - p_mdx) / (2 * dx);
    del_y = (p_pdy - p_mdy) / (2 * dy);
    gradient = [del_x, del_y];
end