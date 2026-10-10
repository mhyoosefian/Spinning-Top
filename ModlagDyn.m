function dz = ModlagDyn(t,z)
dz = z;
global m h Ixx Iyy Izz g Moment

q1 = z(1); q2 = z(2); q3 = z(3);
u1 = z(4); u2 = z(5); u3 = z(6);
u = [u1;u2;u3];
W = W_Ml(q2);
X = X_Ml;
qdot = W*u + X;
dq1 = qdot(1); dq2 = qdot(2); dq3 = qdot(3);
Qprime_ML = Qpr_Ml(Moment);
Mprime_ML = Mpr_Ml(Ixx,Iyy,Izz,h,m); 
Bprime_ML = Bpr_Ml(Ixx,Izz,dq1,dq2,dq3,g,h,m,q2);

A = W*u + X;
B = Mprime_ML\(Qprime_ML - Bprime_ML);


dz = [A;B];