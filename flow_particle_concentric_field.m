%% Flowing Particles: Reproduce Animations in the Supplementary Material
%% IMPORTANT NOTICE TO AVOID CONFUSION: 
% To properly name variables, we name "psi" as the MODULATED POTENTIAL, and
% we name "phi" as the USER-Defined potential

%% define discretization space
[xmin, xmax, ymin, ymax] = deal(0, 1, 0, 1);
corner = [xmin ymin; xmin ymax; xmax ymax; xmax ymin];
xx = corner(:,1);
yy = corner(:,2);
roi = polyshape(xx, yy); % region of interest
tol = 1e-6; % tolerance to decide if a particle is inside the roi

ct_dx = 2.5e-3;
ct_dy = 2.5e-3;
ct_x = xmin : ct_dx : xmax - ct_dx;
ct_y = ymin : ct_dy : ymax - ct_dy;
ct_num_x = length(ct_x);
ct_num_y = length(ct_y);
[ct_X, ct_Y] = meshgrid(ct_x,ct_y);
fprintf("space discretized\n")

%% load contour matrix of psi' for plotting and animation purpose
% concentric_contour_matirx is saved from "contour_prep_concentric_field.m"
load contour/concentric_contour_matrix
fprintf("pre-calculated contour matrix loaded\n")

%% Define obstacles, SC map, psi (psi' in paper), phi (psi in paper), and ramp function 
% Obstacle 1: a circle
[x1, y1] = deal(0.57,0.63);
R1 = 0.1;
theta = 0 : 2 * pi / 200 : 2 * pi;
[vx1, vy1] = deal(R1 * cos(theta) + x1, R1 * sin(theta) + y1);
circle = @(x,y) (x - x1).^2 + (y - y1).^2 - R1^2;

% Obstacle 2: a cross
R2 = 0;
vx2 = [0.55, 0.6, 0.55, 0.5, 0.45, 0.4, 0.45, 0.4, 0.45, 0.5, 0.55, 0.6] + 0.1;
vy2 = [0.5, 0.55, 0.6, 0.55, 0.6, 0.55, 0.5, 0.45, 0.4, 0.45, 0.4, 0.45] - 0.15;
[x2, y2] = deal(mean(vx2), mean(vy2));
proxy2 = boundary_proxy(vx2, vy2);
p2 = polygon(complex(proxy2(:,1), proxy2(:,2)));
f2 = extermap(p2);
g2 = inv(f2);
pgon2 = polyshape(vx2', vy2');

% Obstacle 3: a polygon
edge_num = 8;
d_theta = 2 * pi / edge_num;
theta = 0 : d_theta : 2 * pi - d_theta;
[x3, y3] = deal(0.3, 0.5);
R3 = 0.15;
vx3 = zeros(1, edge_num);
vy3 = zeros(1, edge_num);
for i = 1 : length(theta)
    if rem(i, 2) == 1 
        vx3(i) = R3 * cos(theta(i))+ x3;
        vy3(i) = R3 * sin(theta(i))+ y3;
    else
        vx3(i) = R3 * cos(theta(i)) / 3 + x3;
        vy3(i) = R3 * sin(theta(i)) / 3 + y3;
    end
end
proxy3 = boundary_proxy(vx3, vy3);
p3 = polygon(complex(proxy3(:,1), proxy3(:,2)));
f3 = extermap(p3);
g3 = inv(f3);
pgon3 = polyshape(vx3', vy3');

% Obstacle info:
obj_list = ["circle", "polygon", "polygon"];
obj_num = length(obj_list);
cx = [x1, x2, x3]; % x coordinates of center of polygon
cy = [y1, y2, y3]; % y coordinates of center of polygon
center = {cx, cy};
vx = {"None", [vx2], [vx3]};
vy = {"None", [vy2], [vy3]};
vertex = {vx, vy};
cr = [R1, R2, R3]; % radius of regular obstacles
% SC map info
% consturct the forward map f: interior of disk -> exteriour of polygon
f = {"None", f2, f3};
% construct g: exteriour of polygon -> interior of disk, as the inverse of f.
g = {"None", g2, g3};

% define ramp function
ramp = @(r) 1 .* (r >= 1) + ...
    (15 * r / 8 - 10 * r.^3 / 8 + 3 * r.^5 / 8) .* ((r < 1) & (r > -1)) + ...
    (-1) * (r <= -1);

% define the user-defined potential
phi = @(x,y) sqrt((x - 0.5).^2 + (y - 0.5).^2);

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

%% Initialize particles in random position OR load particles' position
init = 0;
if init
    concentric_pos = seed_particle(obj_list, center, vertex, cr);
    old_pos = concentric_pos;
    save('position/concentric_pos;.mat','concentric_pos;');
    fprintf("particle position seeded\n");  
else
    load position/concentric_pos
    old_pos = concentric_pos;
    fprintf("particle position loaded\n");  
end

%% make animation
d0 = 0.1;
dt = 2e-4;
iter_num = 1500;
offset = 1e-2;
points_num = size(old_pos, 1);

inbox = @(left, right, bottom, top, pos) (left <= pos(:,1)) .* (pos(:,1) < right) .*...
                                             (bottom <= pos(:,2)) .* (pos(:,2) < top);

record_avi = 0;
if(record_avi)
    writerObj = VideoWriter('concentric.avi');
    open(writerObj);
end

fig = figure(10); clf(fig);
[pos1, pos2, pos3] = deal(old_pos, old_pos, old_pos);
for iter = 1 : iter_num
    pos3 = pos2; pos2 = pos1;
    fprintf("iteration num = %d\n", iter);
    % update velocity
    velocity = nabla_particle_potential(psi, d0, old_pos, obj_list, center, vertex, cr, f, g);
    velocity = [velocity(:,2), -velocity(:,1)];
    old_pos = old_pos + velocity * dt;
    % Let particle adhere to isocontour, if needed. Uncomment the following line, but this is time-consuming.
    % old_pos = adhere_isocontour(init_p, psi, d0, old_pos, obj_list, center, vertex, cr, f, g);
    old_pos = tweak_particle(old_pos, obj_list, center, vertex, cr, offset);
    pos1 = old_pos;
    
    % plot frame
    fimplicit(circle); hold on; fill(vx1, vy1, [233,185,136]/255); hold on;
    plot(pgon2); hold on;
    plot(pgon3); hold on;
    plot([proxy2(:,1); proxy2(1,1)], [proxy2(:,2); proxy2(1,2)], "r", 'LineWidth',1); hold on;
    plot([proxy3(:,1); proxy3(1,1)], [proxy3(:,2); proxy3(1,2)], "r", 'LineWidth',1); hold on;
    plot([pos1(:,1)'; pos2(:,1)'; pos3(:,1)'], [pos1(:,2)'; pos2(:,2)'; pos3(:,2)'], "black.-", 'MarkerSize',9,'LineWidth',2);
    hold on;
    contour(ct_X, ct_Y, concentric_contour_matrix, 80); hold on;
    axis equal square
    xlim([xmin, xmax]); ylim([ymin, ymax]);
    xticks(xmin : 0.1 : xmax); yticks(ymin : 0.1 : ymax);
    xlabel("X"); ylabel("Y");
    hold off;
    pause(0.01)
    
    drawnow;
    %get the frame
    X = getframe;
    %save the frame for replay
    F(iter) = X;
    %write the frame to disk
    if(record_avi)
        writeVideo(writerObj,X);
    end
end

if(record_avi)
    close(writerObj);
end
