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

ct_dx = 2.5e-3;
ct_dy = 2.5e-3;
ct_x = xmin : ct_dx : xmax - ct_dx;
ct_y = ymin : ct_dy : ymax - ct_dy;
ct_num_x = length(ct_x);
ct_num_y = length(ct_y);
[ct_X, ct_Y] = meshgrid(ct_x,ct_y);
fprintf("space discretized\n")

%% load contour
load contour/fig6_contour_matrix
fprintf("pre-calculated contour matrix loaded\n")

%% Define obstacles, SC map, psi (psi' in paper), phi (psi in paper), and ramp function 
% Obstacle 1: a convex polygon
edge_num = 8;
d_theta = 2 * pi / edge_num;
theta = 0 : d_theta : 2 * pi - d_theta;
[x1, y1] = deal(0.25, 0.48);
R1 = 0.1;
vx1 = R1 * cos(theta)+ x1;
vy1 = R1 * sin(theta)+ y1;
proxy1 = boundary_proxy(vx1, vy1);
p1 = polygon(complex(proxy1(:,1), proxy1(:,2)));
f1 = extermap(p1);
g1 = inv(f1);
pgon1 = polyshape(vx1', vy1');

% Obstacle 2: a circle
[x2, y2] = deal(0.68, 0.48);
R2 = 0.08;
theta = 0 : 2 * pi / 200 : 2 * pi;
[vx2, vy2] = deal(R2 * cos(theta) + x2, R2 * sin(theta) + y2);
circle2 = @(x,y) (x - x2).^2 + (y - y2).^2 - R2^2;

% Obstacle 3: a polygon boat
R3 = 0;
vx3 = [0.75, 0.60, 0.55, 0.48, 0.46, 0.36, 0.46, 0.64] - 0.05;
vy3 = [0.50, 0.50, 0.62, 0.62, 0.50, 0.50, 0.38, 0.38] - 0.15;
[x3, y3] = deal(mean(vx3), mean(vy3));
proxy3 = boundary_proxy(vx3, vy3);
p3 = polygon(complex(proxy3(:,1), proxy3(:,2)));
f3 = extermap(p3);
g3 = inv(f3);
pgon3 = polyshape(vx3', vy3');

% Obstacle 4: fox face
R4 = 0;
vx4 = [0.67, 0.63, 0.61, 0.59, 0.41, 0.39, 0.37, 0.33, 0.33, 0.37, 0.41, 0.5, 0.59, 0.63, 0.67];
vy4 = [0.56, 0.6, 0.64, 0.6, 0.6, 0.64, 0.6, 0.56, 0.46, 0.42, 0.42, 0.38, 0.42, 0.42, 0.46] + 0.17;
[x4, y4] = deal(mean(vx4), mean(vy4));
proxy4 = boundary_proxy(vx4, vy4);
p4 = polygon(complex(proxy4(:,1), proxy4(:,2)));
f4 = extermap(p4);
g4 = inv(f4);
pgon4 = polyshape(vx4', vy4');

% Obstacle info:
obj_list = ["polygon", "circle", "polygon", "polygon"];
obj_num = length(obj_list);
cx = [x1, x2, x3, x4]; % x coordinates of center of polygon
cy = [y1, y2, y3, y4]; % y coordinates of center of polygon
center = {cx, cy};
vx = {[vx1], "None", [vx3], [vx4]};
vy = {[vy1], "None", [vy3], [vy4]};
vertex = {vx, vy};
cr = [R1, R2, R3, R4]; % radius of regular obstacles
% SC map info
% consturct the forward map f: interior of disk -> exteriour of polygon
f = {f1, "None", f3, f4};
% construct g: exteriour of polygon -> interior of disk, as the inverse of f.
g = {g1, "None", g3, g4};

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

%% Initialize particles in random position OR load particles' position
init = 0;
if init
    fig6_pos = seed_particle(obj_list, center, vertex, cr);
    old_pos = fig6_pos;
    save('position/fig6_pos.mat','fig6_pos');
    fprintf("particle position seeded\n");  
else
    load position/fig6_pos
    old_pos = fig6_pos;
    fprintf("particle position loaded\n");
end

%% make animation
d0 = 0.1;
dt = 2e-4;
iter_num = 1500;
offset = 5e-3;
points_num = size(old_pos, 1);
inbox = @(left, right, bottom, top, pos) (left <= pos(:,1)) .* (pos(:,1) < right) .*...
                                             (bottom <= pos(:,2)) .* (pos(:,2) < top);

record_avi = 0;
if(record_avi)
    writerObj = VideoWriter('fig6_animation.avi');
    open(writerObj);
end

fig = figure(10); clf(fig);
[pos1, pos2, pos3] = deal(old_pos, old_pos, old_pos);
for iter = 1 : iter_num
    pos3 = pos2; pos2 = pos1;
    fprintf("iteration num = %d\n", iter);
    velocity = nabla_particle_potential(psi, d0, old_pos, obj_list, center, vertex, cr, f, g);
    velocity = [velocity(:,2), -velocity(:,1)];
    old_pos = old_pos + velocity * dt;
    % Let particle adhere to isocontour, if needed. Uncomment the following line, but this is time-consuming.
    % old_pos = adhere_isocontour(init_p, psi, d0, old_pos, obj_list, center, vertex, cr, f, g);
    old_pos = tweak_particle(old_pos, obj_list, center, vertex, cr, offset);
    pos1 = old_pos;
    
    % plot frame
    fimplicit(circle2); hold on; fill(vx2, vy2, [233,185,136]/255); hold on;
    plot(pgon1); hold on;
    plot(pgon3); hold on;
    plot(pgon4); hold on;
    plot([proxy1(:,1); proxy1(1,1)], [proxy1(:,2); proxy1(1,2)], "r", 'LineWidth',1.5); hold on;
    plot([proxy3(:,1); proxy3(1,1)], [proxy3(:,2); proxy3(1,2)], "r", 'LineWidth',1.5); hold on;
    plot([proxy4(:,1); proxy4(1,1)], [proxy4(:,2); proxy4(1,2)], "r", 'LineWidth',1.5); hold on;
    plot([pos1(:,1)'; pos2(:,1)'; pos3(:,1)'], [pos1(:,2)'; pos2(:,2)'; pos3(:,2)'], "black.-", 'MarkerSize',9,'LineWidth',2);
    hold on;
    contour(ct_X, ct_Y, fig6_contour_matrix, 100); hold on; caxis([0 4]);
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

