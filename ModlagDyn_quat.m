function dz = ModlagDyn_quat(t,z)
dz = z;
global m h Moment g Ixx Iyy Izz
e0 = z(1); e1 = z(2); e2 = z(3); e3 = z(4);
mag = sqrt(e0^2+e1^2+e2^2+e3^2);
e0 = e0/mag;e1 = e1/mag;e2 = e2/mag;e3 = e3/mag;
u1_quat = z(5); u2_quat = z(6); u3_quat = z(7);

u_quat = [u1_quat;u2_quat;u3_quat];
W_quat = W_Ml_quat(e0,e1,e2,e3);
X_quat = X_Ml_quat;
Mprime_ML_quat = Mpr_Ml_quat(Ixx,Iyy,Izz,h,m);
Qprime_ML_quat = Qpr_Ml_quat(Moment);
Bprime_ML_quat_itu = Bpr_Ml_quat_itu(Ixx,Iyy,Izz,e0,e1,e2,e3,g,h,m,u1_quat,u2_quat,u3_quat);

dz = [W_quat*u_quat+X_quat; Mprime_ML_quat\(Qprime_ML_quat-Bprime_ML_quat_itu)];