%% Local stability scan for 2-DOF planar arm (CoM model: lc1, lc2)
% - Uses CoM-based inertia (lc1, lc2)
% - Optionally uses link rotational inertias I1,I2 (set to 0 to ignore)
% - Computes A,B by finite-difference of the full nonlinear f(x,u)
% - Finds poses where max real(eig(A)) < 0 (locally asymptotically stable)

clear; clc; format long;

% ---------- Parameters (edit as needed) ----------
params.m1  = 1.37;       % kg
params.m2  = 1.16;       % kg
params.l1  = 0.30;      % m
params.l2  = 0.25;      % m
params.lc1 = params.l1/2;   % m, CoM location link1 (default midpoint)
params.lc2 = params.l2/2;   % m, CoM location link2 (default midpoint)
params.I1  = 0.0121;       % kg*m^2 (link rotational inertia about its CoM or joint; set 0 if ignoring)
params.I2  = 0.0074;       % kg*m^2
params.g   = 9.81;      % m/s^2

% ---------- Grid settings ----------
nsteps = 81;                          % grid resolution (81x81)
t1grid = linspace(-pi, pi, nsteps);   % theta1 grid
t2grid = linspace(-pi, pi, nsteps);   % theta2 grid

MaxRealEig = nan(nsteps,nsteps);
CtlRank = nan(nsteps,nsteps);

% ---------- Sweep grid ----------
fprintf('Scanning %dx%d grid of poses ... this may take a while.\n', nsteps, nsteps);
for i=1:nsteps
    for j=1:nsteps
        th1 = t1grid(i);
        th2 = t2grid(j);
        [A_num, B_num] = compute_AB(th1, th2, params);
        ev = eig(A_num);
        MaxRealEig(i,j) = max(real(ev));
        CtrbRank = rank(ctrb(A_num,B_num));
        CtlRank(i,j) = CtrbRank;
    end
end
fprintf('Grid scan complete.\n');

% ---------- Plot heatmap ----------
figure('Name','Max Real(eig(A)) over theta grid','NumberTitle','off');
imagesc(t2grid, t1grid, MaxRealEig); set(gca,'YDir','normal');
xlabel('\theta_2 (rad)'); ylabel('\theta_1 (rad)');
title('Max Real(eig(A)) over theta grid (negative = stable)');
colorbar;
hold on;
% overlay contour where MaxRealEig = 0
contour(t2grid, t1grid, MaxRealEig, [0 0], 'k', 'LineWidth', 1.2);
hold off;

% ---------- Report candidate stable points ----------
[idx_i, idx_j] = find(MaxRealEig < 0);
fprintf('Found %d grid points with all eigenvalues negative (locally stable).\n', numel(idx_i));

% ---------- Most-stable grid pose ----------
[minVal, minIndexLinear] = min(MaxRealEig(:));
[i_min, j_min] = ind2sub(size(MaxRealEig), minIndexLinear);
best_th1 = t1grid(i_min); best_th2 = t2grid(j_min);
fprintf('Most-stable grid pose: theta1=%.6f rad (%.2f deg), theta2=%.6f rad (%.2f deg), maxReEig=%.6e\n', ...
    best_th1, rad2deg(best_th1), best_th2, rad2deg(best_th2), minVal);

% ---------- refine using fminsearch ----------
objfun = @(th) max(real(eig(compute_AB(th(1), th(2), params))));
options = optimset('Display','iter','TolX',1e-8,'TolFun',1e-8);
[xopt, fval, exitflag, output] = fminsearch(objfun, [best_th1, best_th2], options);
fprintf('Refined pose: theta1=%.6f rad (%.2f deg), theta2=%.6f rad (%.2f deg), maxReEig=%.6e\n', ...
    xopt(1), rad2deg(xopt(1)), xopt(2), rad2deg(xopt(2)), fval);

% ---------- check controllability and eigenvalues at refined pose ----------
[Aopt, Bopt] = compute_AB(xopt(1), xopt(2), params);
rankCtrb_opt = rank(ctrb(Aopt, Bopt));
ev_opt = eig(Aopt);
fprintf('At refined pose: controllability rank = %d (n=4)\n', rankCtrb_opt);
fprintf('Eigenvalues at refined pose:\n'); disp(ev_opt);


%% -------------------------
% compute_AB: returns numeric A,B for given theta1,theta2 (thdot=0)
% --------------------------
function [A_num, B_num] = compute_AB(theta1, theta2, params)
    % Unpack params
    m1  = params.m1;  m2  = params.m2;
    l1  = params.l1;  l2  = params.l2;
    lc1 = params.lc1; lc2 = params.lc2;
    I1  = params.I1;  I2  = params.I2;
    g   = params.g;

    % Precompute coupling term
    a = m2 * l1 * lc2;  % standard for CoM-based model

    % Define equilibrium torque (gravity balance) at zero velocity
    % Gravity components (numeric)
    G1 = m1*g*lc1*cos(theta1) + m2*g*(l1*cos(theta1) + lc2*cos(theta1+theta2));
    G2 = m2*g*lc2*cos(theta1+theta2);
    tau0 = [G1; G2];

    % Nonlinear dynamics function f(x,u)
    function xdot = f_nl(x,u)
        th1 = x(1); th2 = x(2); th1d = x(3); th2d = x(4);

        % Inertia matrix (CoM-based) including optional link inertias I1,I2
        M11 = m1*lc1^2 + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(th2)) + I1 + I2;
        M12 = m2*(lc2^2 + l1*lc2*cos(th2)) + I2;
        M21 = M12;
        M22 = m2*lc2^2 + I2;
        Mmat = [M11, M12; M21, M22];

        % Coriolis / centrifugal vector (2x1)
        C1 = -2*a*sin(th2)*th1d*th2d - a*sin(th2)*th2d^2;
        C2 =  a*sin(th2)*th1d^2;
        Cvec_local = [C1; C2];

        % Gravity vector
        G1_loc = m1*g*lc1*cos(th1) + m2*g*(l1*cos(th1) + lc2*cos(th1+th2));
        G2_loc = m2*g*lc2*cos(th1+th2);
        G_local = [G1_loc; G2_loc];

        % accelerations
        thdd_local = Mmat \ (u - Cvec_local - G_local);

        xdot = [th1d; th2d; thdd_local(1); thdd_local(2)];
    end

    % Evaluate Jacobians numerically at x0 = [theta1; theta2; 0; 0], u0 = tau0
    x0 = [theta1; theta2; 0; 0];
    u0 = tau0;

    n = 4; m = 2;
    epsv = 1e-6;

    f0 = f_nl(x0, u0);

    A_num = zeros(n);
    B_num = zeros(n, m);

    % finite difference for A
    for k = 1:n
        dx = zeros(n,1); dx(k) = epsv;
        A_num(:,k) = (f_nl(x0 + dx, u0) - f0) / epsv;
    end

    % finite difference for B
    for k = 1:m
        du = zeros(m,1); du(k) = epsv;
        B_num(:,k) = (f_nl(x0, u0 + du) - f0) / epsv;
    end
end
