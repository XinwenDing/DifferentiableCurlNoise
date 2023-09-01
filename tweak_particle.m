function new_location = tweak_particle(curr_pos, obj_list, center, vertex, cr, offset)
% This function tweaks particle whose initial randomized location is inside an object
% obj_list: list of objects.
% center: e.g. center = {cx, cy}; center of the objects, where cx and cy
%         stores the x- and y-coordinates of the center of the objects
% vertex: e.g. vertex = {vx, vy}; vertices of the objects, where vx and vy
%         stores the x- and y-coordinates of the vertices of the objects
% cr: Radius of objects, if properly defined (e.g. circles and regular polygons have radius)
% offset: if particle is inside any objects, then how far away to drag it out
% For definition and data structure of obj_list, center, vertex, cr, please
% search (Ctrl-F) and check the "Obstacle info" part in any file starts with "flow particle"
    new_location = curr_pos;
    particle_num = size(curr_pos,1);
    obj_num = length(obj_list);
    cx = center{1}; % eg. cx = [x1, x2];
    cy = center{2}; % eg. cy = [y1, y2];
    vx = vertex{1}; % eg. vx = {"None", [vx2]};
    vy = vertex{2}; % eg. vy = {"None", [vy2]};
    for i = 1 : particle_num
        pos = curr_pos(i,:);
        for idx = 1 : obj_num
            if obj_list(idx) == "circle"
                [d, cp] = min_circ_dist(cx(idx), cy(idx), cr(idx), pos(1), pos(2));
                if d < 0
                    dir = cp - pos;
                    unit_dir = dir / norm(dir);
                    new_location(i,:) = pos + (abs(d) + offset) * unit_dir;
                    break;
                end
            elseif obj_list(idx) == "polygon"
                [d, cp, flip] = min_poly_dist(vx{idx}, vy{idx}, pos(1), pos(2));
                if d <= 1e-6
                    dir = cp - pos;
                    unit_dir = dir / norm(dir);
                    if flip
                        new_location(i,:) = pos + (abs(d) + 1e-2) * unit_dir;
                        break;
                    else
                        new_location(i,:) = pos - (abs(d) + offset) * unit_dir;
                        break;
                    end
                end
            end
        end
    end
end