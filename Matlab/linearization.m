% Linearize 2-DOF planar arm with CoM (lc1,lc2) and zero joint inertias

clear; clc; format long;

%% Symbolic variables
syms th1 th2 th1d th2d real
syms m1 m2 l1 l2 lc1 lc2 I1 I2 g real
syms tau1 tau2 real

% state & input vectors
x_sym = [th1; th2; th1d; th2d];
tau_sym = [tau1; tau2];

%% Dynamics using CoM distances (lc1, lc2) and no motor inertias
a_sym = m2*l1*lc2;   % coupling coefficient (note lc2 used)

% Inertia matrix (CoM-based)

% then use these expressions for M
M11 = m1*lc1^2 + m2*(l1^2 + lc2^2 + 2*l1*lc2*cos(th2)) + I1 + I2;
M12 = m2*(lc2^2 + l1*lc2*cos(th2)) + I2;
M21 = M12;
M22 = m2*lc2^2 + I2;
M_sym = [M11, M12; M21, M22];


% Coriolis/centrifugal vector (C*thetadot)
C1 = -2*a_sym*sin(th2)*th1d*th2d - a_sym*sin(th2)*th2d^2;
C2 =  a_sym*sin(th2)*th1d^2;
Cvec_sym = [C1; C2];

% Gravity vector using CoM distances
G1 = m1*g*lc1*cos(th1) + m2*g*(l1*cos(th1) + lc2*cos(th1+th2));
G2 = m2*g*lc2*cos(th1+th2);
G_sym = [G1; G2];

% state derivatives
Minv_sym = simplify(inv(M_sym));
thdd_sym = simplify(Minv_sym*(tau_sym - Cvec_sym - G_sym));
f_sym = [th1d; th2d; thdd_sym(1); thdd_sym(2)];   % f(x,u)

%% Jacobians (symbolic)
A_sym = jacobian(f_sym, x_sym);
B_sym = jacobian(f_sym, tau_sym);

% Outputs: choose mu = y = [theta1; theta2] (mu = y model)
y_sym = [th1; th2];
C_sym = jacobian(y_sym, x_sym);
D_sym = jacobian(y_sym, tau_sym);  % will be zero
%%
% numeric parameters (add I1 and I2)
numvals_pairs = {
    m1, 1.37;
    m2, 1.16;
    l1, 0.30;
    l2, 0.25;
    lc1, 0.15;
    lc2, 0.125;
    I1,  0.0121;    % link inertia (kg*m^2)
    I2,  0.0074;   % link inertia (kg*m^2)
    g, 9.81
};


% Equilibrium pose (radians) and zero velocities
% Example: theta1 = -180 deg -> -pi rad, theta2 = 8.9 deg -> 0.155 rad
th0 = [-3.141593; 0.549779];    % set equilibrium 
thd0 = [0; 0];

% Compute equilibrium torque tau0 = G(th0)
% substitute numeric params into G_sym
subs_list = [th1, th2, numvals_pairs(:,1)'];
subs_vals = [th0(1), th0(2), cell2mat(numvals_pairs(:,2)')];

G0_sym = subs(G_sym, subs_list, subs_vals);
G0 = double(G0_sym);
tau0 = G0(:);

%% Build substitution mapping for A,B,C,D evaluation
% create pairs including state and tau
subs_pairs_full = {
    th1, th0(1);
    th2, th0(2);
    th1d, thd0(1);
    th2d, thd0(2)
};
% append numeric parameter pairs
subs_pairs_full = [subs_pairs_full; numvals_pairs];
% append tau values
subs_pairs_full = [subs_pairs_full; {tau1, tau0(1)}; {tau2, tau0(2)}];

% Convert to vectors for subs
sym_vars = subs_pairs_full(:,1);
sym_vals = subs_pairs_full(:,2);

%% Evaluate A,B,C,D numerically 
A_sub = subs(A_sym, sym_vars, sym_vals);
B_sub = subs(B_sym, sym_vars, sym_vals);
C_sub = subs(C_sym, sym_vars, sym_vals);
D_sub = subs(D_sym, sym_vars, sym_vals);

% Check for any remaining symbolic variables
remA = symvar(A_sub);
remB = symvar(B_sub);
if ~isempty(remA) || ~isempty(remB)
    warning('Some symbolic variables remain:');
    if ~isempty(remA), disp('A has:'); disp(remA); end
    if ~isempty(remB), disp('B has:'); disp(remB); end
end

% Convert to numeric
A_num = double(A_sub);
B_num = double(B_sub);
C_num = double(C_sub);
D_num = double(D_sub);

%% Sanity checks
% residual f(x0,u0) should be near zero
f_sub = subs(f_sym, sym_vars, sym_vals);
f_num = double(f_sub);
fprintf('Residual f(x0,u0) (should be ~0):\n'); disp(f_num);

% eigenvalues
eigA = eig(A_num);
fprintf('Eigenvalues of A:\n'); disp(eigA);
fprintf('Max real(eig(A)) = %+g\n', max(real(eigA)));

% controllability
ctrb_rank = rank(ctrb(A_num, B_num));
fprintf('Controllability rank = %d (need 4 for full control)\n', ctrb_rank);

% Build state-space LTI object
A = A_num; B = B_num; C = C_num; D = D_num;
sys = ss(A,B,C,D);

% Display final matrices
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);
disp('D ='); disp(D);
%%
%Transfer function matrix
G = tf(sys);
G_tf=minreal(G,1e-6);
disp('Transfer function matrix G(s):');
G_tf
