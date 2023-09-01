function potential = get_particle_potential(psi, d0, curr_pos, obj_list, center, vertex, cr, f, g)
% This function calculates the potential value for each query point whose
% coordinats are stored in curr_pos.
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
    cx = center{1}; % eg. cx = [x1, x2];
    cy = center{2}; % eg. cy = [y1, y2];
    particle_num = size(curr_pos, 1);
    obj_idx_list = zeros(particle_num, 1);
    ramp_list = zeros(particle_num, 1);
    complex_cp_list = zeros(particle_num, 1);
    center_list = zeros(particle_num, 1);
    for i = 1 : particle_num
        ptx = curr_pos(i, 1);
        pty = curr_pos(i, 2);
        [ramp_input, obj_idx, cp] = potential_prep(d0, ptx, pty, obj_list, center, vertex, cr);
        ramp_list(i) = ramp_input;
        obj_idx_list(i) = obj_idx;
        center_list(i) = complex(cx(obj_idx), cy(obj_idx));
        if ~isnan(cp)
            complex_cp_list(i) = complex(cp(1), cp(2));
        end
    end
    % Notice that the SC-toolbox we are using cannot accurately calculate
    % the image of SC exterior map and the image of inverse SC exterior map
    % simultaneously. To avoid this problem, we chop our domain [0,1]x[0,1]
    % into multiple (25) smaller squares whose edge is of size 0.2 (or smaller).
    % In each of these smallers squares, we compute SC exterior map and the
    % inverse SC exterior map for multiple query points at a time.
    epsilon = 1e-4;
    edge_size = 0.2;
    square_num = round(1 / edge_size);
    obj_num = length(obj_list);
    % the points whose closest point is on a polygon
    sift_cp = (complex_cp_list == 0);
    target_pos = curr_pos(sift_cp,:);
    target_idx = obj_idx_list(sift_cp);
    all_pgon_cp_list = zeros(sum(sift_cp), 1);
    inbox = @(left, right, bottom, top, pos) (left <= pos(:,1)) .* (pos(:,1) < right) .*...
                                             (bottom <= pos(:,2)) .* (pos(:,2) < top);
    for k = 1 : obj_num
        if obj_list(k) == "polygon"
            sift_obj = (target_idx == k);
            obj_k_pos = target_pos(sift_obj, :);
            curr_pgon_cp_list = zeros(sum(sift_obj),1);
            for i = 1 : square_num
                for j = 1 : square_num
                    % define current square of size 0.2
                    left = edge_size * (i - 1);
                    right = edge_size * i;
                    bottom = edge_size * (j - 1);
                    top = edge_size * j;
                    % select particles inside the current box
                    in_box = inbox(left, right, bottom, top, obj_k_pos);
                    local_idx = (in_box == 1);
                    active_pos = obj_k_pos(local_idx,:);
                    % convert R^2 coordinate as a complex number
                    input = complex(active_pos(:,1), active_pos(:, 2));
                    % find inverse SC exterior map
                    retval = eval(g{k}, input);
                    % Scale epsilon close to the boundary of unit circle
                    disk_boundary = (1 - epsilon) * retval ./ abs(retval);
                    % find SC exterior map
                    curr_pgon_cp_list(local_idx) = eval(f{k}, disk_boundary);
                end
            end
            all_pgon_cp_list(sift_obj) = curr_pgon_cp_list;
        end
    end
    complex_cp_list(sift_cp) = all_pgon_cp_list;
    x_input = curr_pos(:, 1);
    y_input = curr_pos(:, 2);
    potential = psi(x_input, y_input, ramp_list, real(center_list), imag(center_list), ...
                    real(complex_cp_list), imag(complex_cp_list));
end