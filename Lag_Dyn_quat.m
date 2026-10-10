function dz = Lag_Dyn_quat(t,z)
dz = z;
global m h Moment g Ixx Iyy Izz
e0 = z(1); e1 = z(2); e2 = z(3); e3 = z(4);
% Normalizing the quaternions
mag = sqrt(e0^2+e1^2+e2^2+e3^2);
e0 = e0/mag;e1 = e1/mag;e2 = e2/mag;e3 = e3/mag;

de0 = z(5); de1 = z(6); de2 = z(7); de3 = z(8);
dq_quat = [de0;de1;de2;de3];

M_quat = MassMat_quat(Ixx,Iyy,Izz,e0,e1,e2,e3,h,m);
B_quat = Bias_quat(Ixx,Iyy,Izz,de0,de1,de2,de3,e0,e1,e2,e3,g,h,m);
adot_quat = ConsJacobdot_quat(de0,de1,de2,de3);
a_quat = ConsJacob_quat(e0,e1,e2,e3);
Q_quat = GenForc_quat(Moment,e0,e1,e2,e3);

% We use Augmented Method to solve the dynamics of the system
F = Q_quat - B_quat;
A_Au = [M_quat -a_quat.'; -a_quat 0];
B_Au = [F; adot_quat*dq_quat];

dz = [dq_quat; A_Au\B_Au];