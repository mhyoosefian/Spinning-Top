clc; clear all; close all; %#ok<*CLALL>


%% Euler Angles Based Solution
syms q1 q2 q3 real %psi, theta, phi
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms h m g Moment Ixx Iyy Izz real
syms t 
syms u1 u2 u3 real
syms du1 du2 du3 real
syms u1_quat u2_quat u3_quat real
syms du1_quat du2_quat du3_quat real

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];
ddq = [ddq1;ddq2;ddq3];

u = [u1;u2;u3];
du = [du1;du2;du3];

u_quat = [u1_quat;u2_quat;u3_quat];
du_quat = [du1_quat;du2_quat;du3_quat];

Mz = Moment*[0; 0; 1];     %In x2y2z2 coordinate system

%Kinematics
w = [-dq1*sin(q2); dq2; dq1*cos(q2)+dq3];     %In x2y2z2 coordinate system
r_G_O = [0;0;h];     %In x2y2z2 coordinate system
V_G = cross(w, r_G_O);     %In x2y2z2 coordinate system

% Kentic Energy
I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
T = 0.5*m*dot(V_G, V_G) + 0.5*dot(w,I*w);

% Potential Energy
V = m*g*h*cos(q2);

%Lagrangian
L = T - V;

%Mechanical Energy and its derivative
E = T + V;
matlabFunction(E, 'File', 'Energy');
Edot = zeros(size(E));
for i=1:length(q)
    Edot = Edot + diff(E,q(i))*dq(i) + diff(E,dq(i))*ddq(i);
end
Edot = simplify(Edot);
matlabFunction(Edot, 'File', 'Energydot');

%Power
P = dot(Mz, w);
% matlabFunction(P, 'File', 'Power');

% Matrices Calculation
pl_pdq = jacobian(L,dq).';
M = jacobian(pl_pdq,dq);
M = simplify(M);
matlabFunction(M, 'File', 'MassMat');

pl_pq = jacobian(L,q).';
B = diff(pl_pdq,t) + jacobian(pl_pdq,q)*dq - pl_pq;
B = simplify(B);
matlabFunction(B, 'File', 'Bias');

Q1 = diff(w,dq1).'*Mz;
Q2 = diff(w,dq2).'*Mz;
Q3 = diff(w,dq3).'*Mz;
Q = [Q1; Q2; Q3];
% matlabFunction(Q, 'File', 'GenForc');

% Quasi velocities
U = sym('U',[3,1]);
U(1) = w(1);
U(2) = w(2);
U(3) = w(3);
Y = jacobian(U, dq);
Z = U - Y*dq;

% Quasi-related Matrices
W = inv(Y);
X = -Y\Z;

% Constructing qdot in terms of u
dq_itu = W*u+X; %#ok<*MINV>

% Constructing Velocities in terms of u
w_itu = subs(w,{dq(1),dq(2),dq(3)},{dq_itu(1),dq_itu(2),dq_itu(3)});
V_G_itu = subs(V_G,{dq(1),dq(2),dq(3)},{dq_itu(1),dq_itu(2),dq_itu(3)});

% Partial Velocities
w1 = diff(w_itu, u1);
w2 = diff(w_itu, u2);
w3 = diff(w_itu, u3);

V1 = diff(V_G_itu, u1);
V2 = diff(V_G_itu, u2);
V3 = diff(V_G_itu, u3);

% Active Forces
F = -m*g*[-sin(q2);0;cos(q2)];

% Generalized Force vector in quasi-coordinate space
Qprime1 = dot(F,V1) + dot(Mz,w1);
Qprime2 = dot(F,V2) + dot(Mz,w2);
Qprime3 = dot(F,V3) + dot(Mz,w3);
Qprime = [Qprime1;Qprime2;Qprime3];

% Generalized force vector in generalized coordinate space (including
% gravity force)
% QQ1 = Q1 + diff(V_G,dq1).'*F;
% QQ2 = Q2 + diff(V_G,dq2).'*F;
% QQ3 = Q3 + diff(V_G,dq3).'*F;
% QQ = [QQ1;QQ2;QQ3];

%Constructing Generalized Force Vector in quasi-coordinate space using that
%of generalized coordinate space
% Qprime2 = W.'*QQ;

% Modified Lagrange Method
Wdot = zeros(size(W));
Xdot = zeros(size(X));
for i = 1:length(q)
    Wdot = Wdot + diff(W, q(i))*dq(i);
    Xdot = Xdot + diff(X, q(i))*dq(i);
end
Wdot = simplify(Wdot);
Xdot = simplify(Xdot);
Mprime_ML = W.'*M*W;
Bprime_ML = W.'*M*(Wdot*U + Xdot) + W.'*B;
Qprime_ML = W.'*Q;

Mprime_ML = simplify(Mprime_ML);
Bprime_ML = simplify(Bprime_ML);
Qprime_ML = simplify(Qprime_ML);

matlabFunction(Mprime_ML, 'File', 'Mpr_Ml');
matlabFunction(Bprime_ML, 'File', 'Bpr_Ml');
matlabFunction(Qprime_ML, 'File', 'Qpr_Ml');
matlabFunction(W, 'File', 'W_Ml');
matlabFunction(X, 'File', 'X_Ml');
matlabFunction(U, 'File', 'u_Ml');
matlabFunction(Wdot, 'File', 'Wdot_Ml');
matlabFunction(Xdot, 'File', 'Xdot_Ml');

%Constructing ddq in terms of u
Wdot_itu = subs(Wdot, {dq(1),dq(2),dq(3)}, {dq_itu(1),dq_itu(2),dq_itu(3)});
ddq_itu = Wdot_itu*u + W*du + Xdot;
Bprime_ML_itu = subs(Bprime_ML, {dq(1),dq(2),dq(3)}, {dq_itu(1),dq_itu(2),dq_itu(3)});
matlabFunction(Bprime_ML_itu, 'File', 'Bpr_Ml_itu');

% Constructing Energy in terms of u
E_itu = subs(E, {dq(1),dq(2),dq(3)}, {dq_itu(1),dq_itu(2),dq_itu(3)});
Edot_itu = subs(Edot, {dq(1),dq(2),dq(3)}, {dq_itu(1),dq_itu(2),dq_itu(3)});
Edot_itu = subs(Edot_itu, {ddq(1),ddq(2),ddq(3)}, {ddq_itu(1),ddq_itu(2),ddq_itu(3)});
E_itu = simplify(E_itu);
Edot_itu = simplify(Edot_itu);
matlabFunction(E_itu, 'File', 'Energy_itu');
matlabFunction(Edot_itu, 'File', 'Energydot_itu');



%% Quaternions Based Solution
syms e0 e1 e2 e3 real
syms de0 de1 de2 de3 real
syms dde0 dde1 dde2 dde3 real

q_quat = [e0;e1;e2;e3];
dq_quat = [de0;de1;de2;de3];
ddq_quat = [dde0;dde1;dde2;dde3];

% Rotation Matrix
R = [e0^2+e1^2-e2^2-e3^2, 2*(e1*e2+e0*e3), 2*(e1*e3-e0*e2); 2*(e1*e2-e0*e3), e0^2+e2^2-e1^2-e3^2, 2*(e2*e3+e0*e1); 2*(e1*e3+e0*e2), 2*(e2*e3-e0*e1),e0^2+e3^2-e1^2-e2^2];

% Angular Velocity
E_quat = [-e1 e0 e3 -e2; -e2 -e3 e0 e1; -e3 e2 -e1 e0];
w_quat = 2*E_quat*dq_quat;

% Kinetic Energy
V_G_quat = cross(r_G_O, w_quat);
T_quat = 0.5*m*dot(V_G_quat, V_G_quat) + 0.5*dot(w_quat,I*w_quat);

% Potential Energy
V_quat = dot(m*g*R.'*r_G_O,[0;0;1]);

%Lagrangian
L_quat = T_quat - V_quat;

%Mechanical Energy and its derivative
En_quat = T_quat + V_quat;
% matlabFunction(En_quat, 'File', 'Energy_quat');
Endot_quat = zeros(size(En_quat));
for i=1:length(q_quat)
    Endot_quat = Endot_quat + diff(En_quat,q_quat(i))*dq_quat(i) + diff(En_quat,dq_quat(i))*ddq_quat(i);
end
Endot_quat = simplify(Endot_quat);
% matlabFunction(Endot_quat, 'File', 'Energydot_quat');

%Power
P_quat = dot(Mz, w_quat);
% matlabFunction(P_quat, 'File', 'Power_quat');

% Matrices Calculation
pl_pdq_quat = jacobian(L_quat,dq_quat).';
M_quat = jacobian(pl_pdq_quat,dq_quat);
M_quat = simplify(M_quat);
% matlabFunction(M_quat, 'File', 'MassMat_quat');

pl_pq_quat = jacobian(L_quat,q_quat).';
B_quat = diff(pl_pdq_quat,t) + jacobian(pl_pdq_quat,q_quat)*dq_quat - pl_pq_quat;
B_quat = simplify(B_quat);
% matlabFunction(B_quat, 'File', 'Bias_quat');

Q1_quat = diff(w_quat,de0).'*Mz;
Q2_quat = diff(w_quat,de1).'*Mz;
Q3_quat = diff(w_quat,de2).'*Mz;
Q4_quat = diff(w_quat,de3).'*Mz;
Q_quat = [Q1_quat; Q2_quat; Q3_quat; Q4_quat];
% matlabFunction(Q_quat, 'File', 'GenForc_quat');

% Constraint
a_quat = q_quat.';
adot_quat = dq_quat.';
% matlabFunction(a_quat, 'File', 'ConsJacob_quat');
% matlabFunction(adot_quat, 'File', 'ConsJacobdot_quat');

% Quasi velocities
U_quat = sym('U_quat',[3,1]);
U_quat(1) = w_quat(1);
U_quat(2) = w_quat(2);
U_quat(3) = w_quat(3);
Y_quat = jacobian(U_quat, dq_quat);
Z_quat = U_quat - Y_quat*dq_quat;

% Quasi-related Matrices
W_quat = inv([Y_quat; a_quat])*[eye(3,3); zeros(1,3)];
X_quat = -inv([Y_quat; a_quat])*[Z_quat; 0];

% Constructing qdot in terms of u
dq_quat_itu = W_quat*u_quat+X_quat; %#ok<*MINV>

% Modified Lagrange Method
Wdot_quat = zeros(size(W_quat));
Xdot_quat = zeros(size(X_quat));
for i = 1:length(q_quat)
    Wdot_quat = Wdot_quat + diff(W_quat, q_quat(i))*dq_quat(i);
    Xdot_quat = Xdot_quat + diff(X_quat, q_quat(i))*dq_quat(i);
end
Wdot_quat = simplify(Wdot_quat);
Xdot_quat = simplify(Xdot_quat);
Mprime_ML_quat = W_quat.'*M_quat*W_quat;
Bprime_ML_quat = W_quat.'*M_quat*(Wdot_quat*U_quat + Xdot_quat) + W_quat.'*B_quat;
Qprime_ML_quat = W_quat.'*Q_quat;

Mprime_ML_quat = simplify(Mprime_ML_quat);
Bprime_ML_quat = simplify(Bprime_ML_quat);
Qprime_ML_quat = simplify(Qprime_ML_quat);

% matlabFunction(Mprime_ML_quat, 'File', 'Mpr_Ml_quat');
% matlabFunction(Bprime_ML_quat, 'File', 'Bpr_Ml_quat');
% matlabFunction(Qprime_ML_quat, 'File', 'Qpr_Ml_quat');
% matlabFunction(W_quat, 'File', 'W_Ml_quat');
% matlabFunction(X_quat, 'File', 'X_Ml_quat');
% matlabFunction(U_quat, 'File', 'u_Ml_quat');
% matlabFunction(Wdot_quat, 'File', 'Wdot_Ml_quat');
% matlabFunction(Xdot_quat, 'File', 'Xdot_Ml_quat');


%Constructing ddq in terms of u
% Wdot_quat_itu = subs(Wdot_quat, {dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)}, {dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
% ddq_quat_itu = Wdot_quat_itu*u_quat + W_quat*du_quat + Xdot_quat;
% Bprime_ML_quat_itu = subs(Bprime_ML_quat, {dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)}, {dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
% matlabFunction(Bprime_ML_quat_itu, 'File', 'Bpr_Ml_quat_itu');

% Constructing Energy and Power in terms of u
% En_quat_itu = subs(En_quat, {dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)}, {dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
% Endot_quat_itu = subs(Endot_quat, {dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)}, {dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
% Endot_quat_itu = subs(Endot_quat_itu, {ddq_quat(1),ddq_quat(2),ddq_quat(3),ddq_quat(4)}, {ddq_quat_itu(1),ddq_quat_itu(2),ddq_quat_itu(3),ddq_quat_itu(4)});
% P_quat_itu = subs(P_quat, {dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)}, {dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
% En_quat_itu = simplify(En_quat_itu);
% Endot_quat_itu = simplify(Endot_quat_itu);
% P_quat_itu = simplify(P_quat_itu);
% matlabFunction(En_quat_itu, 'File', 'Energy_quat_itu');
% matlabFunction(Endot_quat_itu, 'File', 'Energydot_quat_itu');
% matlabFunction(P_quat_itu, 'File', 'Power_quat_itu');




% Constructing Velocities in terms of u
w_quat_itu = subs(w_quat,{dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)},{dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
V_G_quat_itu = subs(V_G_quat,{dq_quat(1),dq_quat(2),dq_quat(3),dq_quat(4)},{dq_quat_itu(1),dq_quat_itu(2),dq_quat_itu(3),dq_quat_itu(4)});
w_quat_itu = simplify(w_quat_itu);
V_G_quat_itu = simplify(V_G_quat_itu);

% Partial Velocities
w1_quat = diff(w_quat_itu, u1_quat);
w2_quat = diff(w_quat_itu, u2_quat);
w3_quat = diff(w_quat_itu, u3_quat);

V1_quat = diff(V_G_quat_itu, u1_quat);
V2_quat = diff(V_G_quat_itu, u2_quat);
V3_quat = diff(V_G_quat_itu, u3_quat);

% Active Forces
F_quat = R*(-m*g*[0;0;1]);

% Generalized Force vector in quasi-coordinate space
Qprime1_quat = dot(F_quat,V1_quat) + dot(Mz,w1_quat);
Qprime2_quat = dot(F_quat,V2_quat) + dot(Mz,w2_quat);
Qprime3_quat = dot(F_quat,V3_quat) + dot(Mz,w3_quat);
Qprime_quat = [Qprime1_quat;Qprime2_quat;Qprime3_quat];
Qprime_quat = simplify(Qprime_quat);

% Generalized force vector in generalized coordinate space (including
% gravity force)
QQ1_quat = Q1_quat + diff(V_G_quat,de0).'*F_quat;
QQ2_quat = Q2_quat + diff(V_G_quat,de1).'*F_quat;
QQ3_quat = Q3_quat + diff(V_G_quat,de2).'*F_quat;
QQ4_quat = Q4_quat + diff(V_G_quat,de3).'*F_quat;
QQ_quat = [QQ1_quat;QQ2_quat;QQ3_quat;QQ4_quat];

%Constructing Generalized Force Vector in quasi-coordinate space using that
%of generalized coordinate space
Qprime2_quat = W_quat.'*QQ_quat;
Qprime2_quat = simplify(Qprime2_quat);