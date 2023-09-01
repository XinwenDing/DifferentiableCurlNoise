function [dist, location] = LSE_distance(vx, vy, x, y)
% This function calculates the LSE distance proposed by Madan and Levin[2022]
% vx and vy are the x- and y-coordinates of vertices of the polygon,
% repectively. x and y are the x- and y-cooridnates of query points.
    alpha = 200;
    location = zeros(1, 2);
    vertex_num = length(vx);
    xq = [x, y];
    distances = zeros(1, vertex_num);
    edge_length = zeros(1, vertex_num);
    closest_locations = zeros(vertex_num, 2);
    lambda = zeros(1, vertex_num);
    if length(vx) == length(vy)
        for i = 1 : vertex_num
            p1 = [vx(i), vy(i)];
            if i == vertex_num
                p2 = [vx(1), vy(1)];
            else
                p2 = [vx(i+1), vy(i+1)];
            end
            % find the closest point from the query point to the line 
            % extended from the current edge
            dp = p2 - p1;
            edge_length(i) = norm(dp);
            xq_p1 = xq - p1;
            if dp(1) == 0
                p0 = [p1(1), xq(2)];
            elseif dp(2) == 0
                p0 = [xq(1), p1(2)];
            elseif norm(cross([dp,0], [xq_p1, 0])) >= 1e-7
                % y - y1 = k(x - x1) ==> y = kx + y1 - kx1
                k1 = dp(2) / dp(1);
                b1 = p1(2) - k1 * p1(1);
                % y = kx + b  ==> b = yq - k * xq, where query = [xq, yq]
                k2 = -1 / k1;
                b2 = y - k2 * x;
                M = [-k1, 1; -k2, 1];
                b = [b1; b2];
                p0 = (M \ b)'; % x0 is on the line of the current edge
            else
                p0 = xq;
            end
            unit_dp = dp / norm(dp);
            s = dot(unit_dp, p2-p0) / norm(dp);
            % locate the closest point somewhere on the edge
            s = max(min(s, 1), 0);
            % the closest point on an edge (including endpoints) to xq 
            x0 = s * p1 + (1 - s) * p2;
            lambda(i) = s;
            distances(i) = norm(xq - x0);
            closest_locations(i,:) = x0;
        end
        round_dist = round(distances, 10);
        [~, idx] = min(round_dist);
        location = closest_locations(idx,:);
        alpha_u = max(1 / min(edge_length), 1200);
        S = alpha / max(alpha, alpha_u);
        A = 2;
        w = @(x) 8 * x.^4 - 16 * x.^3 + 8 * x.^2 + 1/2;
        % Below is a weight function, w, whose second order derivatives at
        % 0 and 1 are 0. Such weight function was discussed in the last
        % section of our supplementary material.
        % w = @(x) -32 * x.^6 + 96 * x.^5 - 96 * x.^4 + 32 * x.^3 + 1/2;
        weight = (A * w(lambda)).^S;
        logsumexp = log(dot(weight, exp(-alpha * distances)));
        dist = - logsumexp / alpha;
    end
end