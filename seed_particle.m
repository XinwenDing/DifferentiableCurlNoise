function init_pos = seed_particle(obj_list, center, vertex, cr)
% This function seeds particle inside the domain of interest, and move them out if the
% initial random genelization of particles locates any of them inside an object
% obj_list: list of objects.
% center: e.g. center = {cx, cy}; center of the objects, where cx and cy
%         stores the x- and y-coordinates of the center of the objects
% vertex: e.g. vertex = {vx, vy}; vertices of the objects, where vx and vy
%         stores the x- and y-coordinates of the vertices of the objects
% cr: Radius of objects, if properly defined (e.g. circles and regular polygons have radius)
% For definition and data structure of obj_list, center, vertex, cr, please
% search (Ctrl-F) and check the "Obstacle info" part in any file starts with "flow particle"
    [xmin, xmax, ymin, ymax] = deal(0, 1, 0, 1);
    [dx, dy] = deal(1e-1, 1e-1);
    x = xmin : dx : xmax;
    y = ymin : dy : ymax;
    [node_num_x, node_num_y]  = deal(length(x), length(y));
    init_pos = zeros(0,0);
    pts_per_subbox = 2;
    for i = 2 : node_num_x
        for j = 2 : node_num_y
            a = [x(i-1), y(j-1)];
            b = [x(i), y(j)];
            pts = (b - a) .* rand(pts_per_subbox, 2) + a;
            init_pos = [init_pos ; pts];
        end
    end
    % tweak random particles outside objects
    [particle_num, ~] = size(init_pos);
    obj_num = length(obj_list);
    cx = center{1}; % eg. cx = [x1, x2];
    cy = center{2}; % eg. cy = [y1, y2];
    vx = vertex{1}; % eg. vx = {"None", [vx2]};
    vy = vertex{2}; % eg. vy = {"None", [vy2]};
    for i = 1 : particle_num
        pos = init_pos(i,:);
        for idx = 1 : obj_num
            if obj_list(idx) == "circle"
                [d, cp] = min_circ_dist(cx(idx), cy(idx), cr(idx), pos(1), pos(2));
                if d < 0
                    dir = cp - pos;
                    unit_dir = dir / norm(dir);
                    init_pos(i,:) = pos + (abs(d) + 2e-2) * unit_dir;
                    break;
                end
            elseif obj_list(idx) == "polygon"
                [d, cp, flip] = min_poly_dist(vx{idx}, vy{idx}, pos(1), pos(2));
                if d <= 1e-6
                    dir = cp - pos;
                    unit_dir = dir / norm(dir);
                    if flip
                        init_pos(i,:) = pos + (abs(d) + 2e-2) * unit_dir;
                        break;
                    else
                        init_pos(i,:) = pos - (abs(d) + 2e-2) * unit_dir;
                        break;
                    end
                end
            end
        end
    end
end