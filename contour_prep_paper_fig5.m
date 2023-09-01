%% contour of psi' for Figure 5 in our paper
% This script calculates the contour matrix of psi' when the user-defined
% potential, psi, is horizontal. The contour matrix is stored in the 
% "contour" foldere
%% IMPORTANT NOTICE TO AVOID CONFUSION: 
% To properly name variables, we name "psi" as the MODULATED POTENTIAL, and
% we name "phi" as the USER-Defined potential

%% define discretization space
[xmin, xmax, ymin, ymax] = deal(0, 1, 0, 1);
corner = [xmin ymin; xmin ymax; xmax ymax; xmax ymin];
xx = corner(:,1);
yy = corner(:,2);
roi = polyshape(xx, yy); % region of interest

ct_dx = 2.5e-3;
ct_dy = 2.5e-3;
ct_x = xmin : ct_dx : xmax - ct_dx;
ct_y = ymin : ct_dy : ymax - ct_dy;
ct_num_x = length(ct_x);
ct_num_y = length(ct_y);
[ct_X, ct_Y] = meshgrid(ct_x,ct_y);
fprintf("space discretized\n")

%% Define obstacles, SC map, psi (psi' in paper), phi (psi in paper), and ramp function 
% Obstacle 1: a convex polygon
edge_num = 6;
d_theta = 2 * pi / edge_num;
theta = 0 : d_theta : 2 * pi - d_theta;
[x1, y1] = deal(0.3, 0.46);
R1 = 0.15;
vx1 = R1 * cos(theta)+ x1;
vy1 = R1 * sin(theta)+ y1;
proxy1 = boundary_proxy(vx1, vy1);
p1 = polygon(complex(proxy1(:,1), proxy1(:,2)));
f1 = extermap(p1);
g1 = inv(f1);
pgon1 = polyshape(vx1', vy1');

% Obstacle 2: a circle
[x2, y2] = deal(0.55, 0.3);
R2 = 0.1;
theta = 0 : 2 * pi / 200 : 2 * pi;
[vx2, vy2] = deal(R2 * cos(theta) + x2, R2 * sin(theta) + y2);
circle2 = @(x,y) (x - x2).^2 + (y - y2).^2 - R2^2;

% Obstacle 3: a plus sign
R3 = 0;
vx3 = [0.62, 0.54, 0.54, 0.45, 0.45, 0.38, 0.38, 0.46, 0.46, 0.54, 0.54, 0.62] + 0.1;
vy3 = [0.54, 0.54, 0.62, 0.62, 0.54, 0.54, 0.46, 0.46, 0.38, 0.38, 0.46, 0.46] + 0.1;
[x3, y3] = deal(mean(vx3), mean(vy3));
proxy3 = boundary_proxy(vx3, vy3);
p3 = polygon(complex(proxy3(:,1), proxy3(:,2)));
f3 = extermap(p3);
g3 = inv(f3);
pgon3 = polyshape(vx3', vy3');

% Obstacle info:
obj_list = ["polygon", "circle",  "polygon"];
obj_num = length(obj_list);
cx = [x1, x2, x3]; % x coordinates of center of polygon
cy = [y1, y2, y3]; % y coordinates of center of polygon
center = {cx, cy};
vx = {[vx1], "None", [vx3]};
vy = {[vy1], "None", [vy3]};
vertex = {vx, vy};
cr = [R1, R2, R3]; % radius of regular obstacles
% SC map info
% consturct the forward map f: interior of disk -> exteriour of polygon
f = {f1, "None", f3};
% construct g: exteriour of polygon -> interior of disk, as the inverse of f.
g = {g1, "None", g3};


% define ramp function
ramp = @(r) 1 .* (r >= 1) + ...
    (15 * r / 8 - 10 * r.^3 / 8 + 3 * r.^5 / 8) .* ((r < 1) & (r > -1)) + ...
    (-1) * (r <= -1);

% define potential
phi = @(x,y) 4 * y;

% define psi' (variable name: psi) to modulate the user-defined potential
% u,v is the coordinate for the closest boundary point
% xc, yc is the coordinate for the center of the obstacle
% x,y is the coordinate for the current point we are interested in
% In our paper, the following psi is the psi' described by Equation 4
% Uncomment the next line when necessary
% psi = @(x, y, r, xc, yc, u, v) ramp(r) * phi(x, y) + (1 - ramp(r)) * phi(xc, yc);
% In our paper, the following psi is the psi' described by Equation 5
% Comment out the next line when necessary
psi = @(x, y, r, xc, yc, u, v) phi(x, y) + (1 - ramp(r)) .* (phi(xc, yc) - phi(u, v));
fprintf("Objects, SC map, phi(psi in paper), psi(psi' in paper), ramp function defined.\n")

%% Calculate the Potential (psi') Contour Matrix
edge_size = 0.2;
square_num = round(xmax / edge_size);
[batch_point_x, batch_point_y] = deal(round(edge_size / ct_dx), round(edge_size / ct_dy));
d0 = 0.1;
epsilon = 1e-4;
fprintf("iteration starts\n");
start_collect = tic;
obj_idx_matrix = zeros(ct_num_x, ct_num_y);
ramp_matrix = zeros(ct_num_x, ct_num_y);
complex_cp_matrix = zeros(ct_num_x, ct_num_y);
center_matrix = zeros(ct_num_x, ct_num_y);
for i = 1 : ct_num_x
    for j = 1 : ct_num_y
        ptx = ct_x(i);
        pty = ct_y(j);
        [ramp_input, obj_idx, cp] = potential_prep(d0, ptx, pty, obj_list, center, vertex, cr);
        ramp_matrix(j,i) = ramp_input;
        obj_idx_matrix(j,i) = obj_idx;
        center_matrix(j,i) = complex(cx(obj_idx), cy(obj_idx));
        if ~isnan(cp)
            complex_cp_matrix(j,i) = complex(cp(1), cp(2));
        end
    end
end
end_collect = toc(start_collect);
fprintf("Time for LSE distance & closest object info collect: %f min\n", end_collect / 60);
fprintf("elements collected!\n");
sift_cp = (complex_cp_matrix == 0);
[target_X, target_Y] = deal(ct_X(sift_cp), ct_Y(sift_cp));
target_idx = obj_idx_matrix(sift_cp);
all_pgon_cp_list = zeros(sum(sum(sift_cp)), 1);
inbox = @(left, right, bottom, top, x, y) (left <= x) .* (x < right) .*...
                                         (bottom <= y) .* (y < top);
start_sc = tic;
for k = 1 : obj_num
    if obj_list(k) == "polygon"
        sift_obj = (target_idx == k);
        [obj_k_X, obj_k_Y] = deal(target_X(sift_obj), target_Y(sift_obj));
        curr_pgon_cp_list = zeros(sum(sum(sift_obj)),1);
        for i = 1 : square_num
            for j = 1 : square_num
                [left, right] = deal(edge_size * (i - 1), edge_size * i);
                [bottom, top] = deal(edge_size * (j - 1), edge_size * j);
                in_box = inbox(left, right, bottom, top, obj_k_X, obj_k_Y);
                local_idx = (in_box == 1);
                [active_pos_X, active_pos_Y] = deal(obj_k_X(local_idx), obj_k_Y(local_idx));
                input = complex(active_pos_X, active_pos_Y);
                retval = eval(g{k}, input);
                disk_boundary = (1 - epsilon) * retval ./ abs(retval);
                curr_pgon_cp_list(local_idx) = eval(f{k}, disk_boundary);
            end
        end
        all_pgon_cp_list(sift_obj) = curr_pgon_cp_list;
    end
    fprintf("complete %dth object\n", k);
end
end_sc = toc(start_sc);
fprintf("Time for SC transform: %f min\n", end_sc / 60);
complex_cp_matrix(sift_cp) = all_pgon_cp_list;
contour_matrix = psi(ct_X, ct_Y, ramp_matrix, real(center_matrix), imag(center_matrix), ...
                     real(complex_cp_matrix), imag(complex_cp_matrix));
fprintf("contour matrix complete\n");

%% plot
figure(1)
fimplicit(circle2); hold on;
plot(pgon1); hold on;
plot(pgon3); hold on;
plot([proxy1(:,1); proxy1(1,1)], [proxy1(:,2); proxy1(1,2)], "r", 'LineWidth',1); hold on;
plot([proxy3(:,1); proxy3(1,1)], [proxy3(:,2); proxy3(1,2)], "r", 'LineWidth',1); hold on;
contour(ct_X, ct_Y, contour_matrix, 80);
title("Curl Noise Potential")
axis equal square
xlim([xmin, xmax]);
ylim([ymin, ymax]);
xticks(xmin : 0.1 : xmax);
yticks(ymin : 0.1 : ymax);
xlabel("X");
ylabel("Y");
hold off;
pause(0.01)
% One can save the variable "contour_matrix" when necessary.
% This matrix is saved as "concentric_contour_matrix" in the "contour" folder.
% Data saving codes are as follows:
% concentric_contour_matrix = contour_matrix;
% save("contour/fig5_contour_matrix.mat", "fig5_contour_matrix");