function dz = Lag_Dyn(t,z)
dz = z;
global m h g Moment Ixx Iyy Izz

q1 = z(1); q2 = z(2); q3 = z(3);
dq1 = z(4); dq2 = z(5); dq3 = z(6);
dq = [dq1;dq2;dq3];

M = MassMat(Ixx,Iyy,Izz,h,m,q2);
B = Bias(Ixx,Izz,dq1,dq2,dq3,g,h,m,q2);
Q = GenForc(Moment,q2);

dz = [dq;M\(Q-B)];