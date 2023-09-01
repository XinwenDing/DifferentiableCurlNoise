function proxy = boundary_proxy(vx, vy)
% This function finds the LSE bondary proxy of polygons whose vertices are
% (vx, vy).
    % first, calculate the outward normal of each edge of the polygon
    normal = zeros(length(vx), 2);
    for i = 1 : length(vx)
        % p1, p2 are the start and the end of the line segment, repectively
        p1 = [vx(i), vy(i)]; 
        if i == length(vx)
            p2 = [vx(1), vy(1)];
        else
            p2 = [vx(i + 1), vy(i + 1)];
        end
        diff = p2 - p1;
        n = diff / norm(diff);
        normal(i,:) = [n(2), -n(1)];
    end
    normal = [normal(end,:); normal];
    pts = @(lambda, start_x, start_y, end_x, end_y)...
                [lambda .* start_x + (1 - lambda) .* end_x, ...
                 lambda .* start_y + (1 - lambda) .* end_y];
    lambda = (1 : -0.05 : 0.05)';
    sample = zeros(0, 0);
    for i = 1 : length(vx)
        % p1, p2 are the start and the end of the line segment, repectively
        p1 = [vx(i), vy(i)]; 
        if i == length(vx)
            p2 = [vx(1), vy(1)];
        else
            p2 = [vx(i + 1), vy(i + 1)];
        end
        points = pts(lambda, p1(1), p1(2), p2(1), p2(2));
        [prev_n, curr_n] = deal(normal(i,:),normal(i+1,:));
        avg_n = prev_n + curr_n;
        avg_n = avg_n / norm(avg_n);
        points(1,:) = points(1,:) + 1e-3 * avg_n;
        points(2:end,:) = points(2:end,:) + 1e-3 * curr_n;
        sample = [sample; points];
    end
    
    %%{
    center = [mean(vx), mean(vy)];
    proxy = zeros(size(sample));
    for i = 1 : length(sample)
        pos = sample(i,:);
        gradient = nabla_LSE_distance(vx, vy, pos(1), pos(2));
        %curr_dist = LSE_distance(vx, vy, pos(1), pos(2));
        [curr_dist, ~] = min_poly_dist(vx, vy, pos(1), pos(2));
        %[curr_dist, gradient, ~] = LSE_dist_grad(vx, vy, pos(1), pos(2));
        %init_dir = pos - center;
        %gradient = init_dir / norm(init_dir);
        norm_square = norm(gradient)^2;
        target = 0;
        diff_dist = curr_dist - target;
        %disp(diff_dist);
        k = 1;
        diff_min = 100;
        best_pos = pos;
        while abs(diff_dist) > 1e-10 && k <= 50
            % improved newton
            pos = pos - diff_dist * gradient / norm_square;
            gradient = nabla_LSE_distance(vx, vy, pos(1), pos(2));
            %curr_dist = LSE_distance(vx, vy, pos(1), pos(2));
            [curr_dist, ~] = min_poly_dist(vx, vy, pos(1), pos(2));
            %[curr_dist, gradient, ~] = LSE_dist_grad(vx, vy, pos(1), pos(2));
            norm_square = norm(gradient)^2;
            diff_dist = curr_dist - target;

            if abs(diff_dist) < diff_min
                best_pos = pos;
                diff_min = abs(diff_dist);
            end
            k = k + 1;
        end
        proxy(i,:) = best_pos;
    end
    %}
    
    %%{
    proxy_size = size(proxy);
    tot_num = proxy_size(1);
    threshold_angle = 175;
    cos_threshold = cos(threshold_angle * pi / 180);
    removed = 0;
    i = 1;
    while i <= tot_num
        if i == 1
            prev_vec = proxy(end, :) - proxy(i,:);
            next_vec = proxy(i+1, :) - proxy(i,:);
        elseif i == tot_num
            prev_vec = proxy(end-1, :) - proxy(end,:);
            next_vec = proxy(1, :) - proxy(i,:);
        else
            prev_vec = proxy(i-1, :) - proxy(i,:);
            next_vec = proxy(i+1, :) - proxy(i,:);
        end
        cos_angle = dot(prev_vec, next_vec) / (norm(prev_vec) * norm(next_vec));
        if cos_angle <= cos_threshold
            proxy(i,:) = [];
            tot_num = tot_num - 1;
            removed = removed + 1;
            %fprintf("i = %d, points removed: %d\n", i, removed);
        else
            i = i + 1;
        end
    end
    %}
end