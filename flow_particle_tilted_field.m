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
tol = 1e-4; % tolerance to decide if a particle is inside the roi

dx = 1e-1;
dy = 1e-1;
x = xmin : dx : xmax;
y = ymin : dy : ymax;
node_num_x = length(x);
node_num_y = length(y);
[X, Y] = meshgrid(x,y);
fprintf("space discretized\n")
%% load contour
load contour/tilted_contour_matrix
fprintf("pre-calculated contour matrix loaded\n")

%% Define obstacles, SC map, psi (psi' in paper), phi (psi in paper), and ramp function 
% Obstacle 1: a polygon crown
R1 = 0;
vx1 = [0.5, 0.36, 0.3, 0.24, 0.1, 0.2, 0.4] + 0.2; %+0.05;
vy1 = [0.36, 0.3, 0.48, 0.3, 0.36, 0.2, 0.2] + 0.05;
[x1, y1] = deal(mean(vx1), mean(vy1));
proxy1 = boundary_proxy(vx1, vy1);
p1 = polygon(complex(proxy1(:,1), proxy1(:,2)));
f1 = extermap(p1);
g1 = inv(f1);
pgon1 = polyshape(vx1', vy1');

% Obstacle 2: a heart
R2 = 0;
vx2 = [0.66, 0.58, 0.5, 0.42, 0.34, 0.34, 0.5, 0.66] - 0.2;
vy2 = [0.5, 0.55, 0.5, 0.55, 0.5, 0.4, 0.3, 0.4] + 0.2;
[x2, y2] = deal(mean(vx2), mean(vy2));
proxy2 = boundary_proxy(vx2, vy2);
p2 = polygon(complex(proxy2(:,1), proxy2(:,2)));
f2 = extermap(p2);
g2 = inv(f2);
pgon2 = polyshape(vx2', vy2');

% Obstacle 3: a polygon
edge_num = 10;
d_theta = 2 * pi / edge_num;
theta = 0 : d_theta : 2 * pi - d_theta;
[x3, y3] = deal(0.7, 0.58);
R3 = 0.12;
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
obj_list = ["polygon", "polygon", "polygon"];
obj_num = length(obj_list);
cx = [x1, x2, x3]; % x coordinates of center of polygon
cy = [y1, y2, y3]; % y coordinates of center of polygon
center = {cx, cy};
vx = {[vx1], [vx2], [vx3]};
vy = {[vy1], [vy2], [vy3]};
vertex = {vx, vy};
cr = [R1, R2, R3]; % radius of regular obstacles
% SC map info
% consturct the forward map f: interior of disk -> exteriour of polygon
f = {f1, f2, f3};
% construct g: exteriour of polygon -> interior of disk, as the inverse of f.
g = {g1, g2, g3};


% define ramp function
ramp = @(r) 1 .* (r >= 1) + ...
    (15 * r / 8 - 10 * r.^3 / 8 + 3 * r.^5 / 8) .* ((r < 1) & (r > -1)) + ...
    (-1) * (r <= -1);

% define potential
phi = @(x,y) x + y;

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
    tilted_pos = seed_particle(obj_list, center, vertex, cr);
    old_pos = tilted_pos;
    save('position/tilted_pos.mat','tilted_pos');
else
    load position/tilted_pos
    old_pos = tilted_pos;
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
    writerObj = VideoWriter('tilted.avi');
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
    plot(pgon1); hold on;
    plot(pgon2); hold on;
    plot(pgon3); hold on;
    contour(ct_X, ct_Y, ex5_contour_matrix, 100); hold on; caxis([0 2]);
    plot([proxy1(:,1); proxy1(1,1)], [proxy1(:,2); proxy1(1,2)], "r", 'LineWidth',1); hold on;
    plot([proxy2(:,1); proxy2(1,1)], [proxy2(:,2); proxy2(1,2)], "r", 'LineWidth',1); hold on;
    plot([proxy3(:,1); proxy3(1,1)], [proxy3(:,2); proxy3(1,2)], "r", 'LineWidth',1); hold on;
    plot([pos1(:,1)'; pos2(:,1)'; pos3(:,1)'], [pos1(:,2)'; pos2(:,2)'; pos3(:,2)'], "black.-", 'MarkerSize',9,'LineWidth',2);
    hold on;
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
