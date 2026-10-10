clc;
clear all; %#ok<*CLALL>
close all;

global m h Moment g Ixx Iyy Izz
m = 2; h =0.3; Moment = 0; g = 9.81; Ixx = 10; Iyy = 10; Izz = 10;

simulation_time = [0:0.01:5]; %#ok<*NBRAK>

%_________________________Choosing Method__________________________________
L = 0;     %1 to use Lagrange Method
ML = 1;     %1 to use Modified Lagrange Method
L_quat = 0;     %1 to use Lagrange Method (Using Quaternions)
ML_quat = 0;     %1 to use Modified Lagrange Method (Using Quaternions)

%% ________________________Euler Based Solutions___________________________


%___________________________Lagrange Method________________________________
if L==1
z0 = [pi/12;pi/4;pi/10;0;0;1];     %States: Z = [q;qdot]
options=odeset('maxstep',0.001);
[t,z] = ode45(@Lag_Dyn,simulation_time,z0, options);
q1 = z(:,1);q2 = z(:,2);q3 = z(:,3);
dq1 = z(:,4);dq2 = z(:,5);dq3 = z(:,6);
figure;
hold on;
plot(t,q1 ,'r-','linewidth',4)
plot(t,q2 ,'b-','linewidth',4)
plot(t,q3 ,'k-','linewidth',4)
legend('q_1','q_2','q_3');
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
ylabel('Generalized Coordinates (rad) ','fontsize',14,'fontweight','bold');

if Moment == 0
E = Energy(Ixx,Iyy,Izz,dq1,dq2,dq3,g,h,m,q2);
figure;
plot(t,100*(E-E(1))/E(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
else
Edot = zeros(length(simulation_time),1);
for i = 1:length(simulation_time)
    M = MassMat(Ixx,Iyy,Izz,h,m,q2(i));
    B = Bias(Ixx,Izz,dq1(i),dq2(i),dq3(i),g,h,m,q2(i));
    Q = GenForc(Moment,q2(i));
    ddq = M\(Q-B);
    Edot(i) = Energydot(Ixx,Iyy,Izz,ddq(1),ddq(2),ddq(3),dq1(i),dq2(i),dq3(i),g,h,m,q2(i)); 
end
P = Power(Moment,dq1,dq3,q2);
figure;
plot(t,(Edot-P),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Power error','fontsize',14,'fontweight','bold');
end
end


%_____________________Modified Lagrange Method_____________________________
if ML == 1
z0 = [pi/12;pi/4;pi/10;0;0;1];     %States: Z = [q;u]
options=odeset('maxstep',0.001);
[t,z] = ode45(@ModlagDyn,simulation_time,z0, options);
q1_Ml = z(:,1); q2_Ml = z(:,2); q3_Ml = z(:,3);
u1 = z(:,4); u2 = z(:,5); u3 = z(:,6);
qdot_Ml = zeros(3,length(simulation_time));
for i=1:length(simulation_time)
u = [u1(i);u2(i);u3(i)];
W = W_Ml(q2_Ml(i));
X = X_Ml;
qdot_Ml(:,i) = W*u + X;
end
figure;
hold on;
plot(t,q1_Ml,'r-','linewidth',4)
plot(t,q2_Ml,'b-','linewidth',4)
plot(t,q3_Ml,'k-','linewidth',4)
legend('\theta_1','\theta_2','\theta_3');
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
ylabel('Generalized Coordinates (rad) ','fontsize',14,'fontweight','bold');

if Moment == 0
E_Ml = Energy(Ixx,Iyy,Izz,qdot_Ml(1,:)',qdot_Ml(2,:)',qdot_Ml(3,:)',g,h,m,q2_Ml);
figure;
plot(t,100*(E_Ml-E_Ml(1))/E_Ml(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
else
Edot_Ml = zeros(length(simulation_time),1);
for i = 1:length(simulation_time)
    u = [u1(i);u2(i);u3(i)];
    Qprime_ML = Qpr_Ml(Moment);
    Mprime_ML = Mpr_Ml(Ixx,Iyy,Izz,h,m); 
    Bprime_ML = Bpr_Ml_itu(Ixx,Izz,g,h,m,q2_Ml(i),u1(i),u2(i),u3(i));
    udot = Mprime_ML\(Qprime_ML-Bprime_ML);
    Edot_Ml(i) = Energydot_itu(Ixx,Iyy,Izz,udot(1),udot(2),udot(3),g,h,m,q2_Ml(i),u1(i),u2(i),u3(i));
end
P_Ml = Power(Moment,qdot_Ml(1,:)',qdot_Ml(3,:)',q2_Ml);
figure;
plot(t,(Edot_Ml-P_Ml),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Power error','fontsize',14,'fontweight','bold');
end
end







%% ____________________Quaternion Based Solutions__________________________



%_________________________Lagrange Method__________________________________
if L_quat==1
z0 = [1;0;0;0;0;2;-1;1;0];     %States: Z = [e0;e1;e2;e3;de0;de1;de2;de3;lamda_quat]
options=odeset('maxstep',0.001);
[t,z] = ode45(@Lag_Dyn_quat,simulation_time,z0, options);
e0 = z(:,1);e1 = z(:,2);e2 = z(:,3);e3 = z(:,4);
de0 = z(:,5);de1 = z(:,6);de2 = z(:,7);de3 = z(:,8);
figure;
hold on;
plot(t,e0 ,'r-','linewidth',4)
plot(t,e1 ,'b-','linewidth',4)
plot(t,e2 ,'k-','linewidth',4)
plot(t,e3 ,'g-','linewidth',4)
legend('e_0','e_1','e_2','e_3');
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
ylabel('Quaternions ','fontsize',14,'fontweight','bold');
%Constraint
f = zeros(length(simulation_time),1);
for i=1:length(simulation_time)
a_quat = ConsJacob_quat(e0(i),e1(i),e2(i),e3(i));
dq_quat = [de0(i);de1(i);de2(i);de3(i)];
f(i,1) = a_quat*dq_quat;
end
figure;
plot(t,f,'k-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Error of Constraint','fontsize',14,'fontweight','bold');

if Moment == 0
En_quat = Energy_quat(Ixx,Iyy,Izz,de0,de1,de2,de3,e0,e1,e2,e3,g,h,m);
figure;
plot(t,100*(En_quat-En_quat(1))/En_quat(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
else
Endot_quat = zeros(length(simulation_time),1);
for i = 1:length(simulation_time)
    dq_quat = [de0(i);de1(i);de2(i);de3(i)];
    M_quat = MassMat_quat(Ixx,Iyy,Izz,e0(i),e1(i),e2(i),e3(i),h,m);
    B_quat = Bias_quat(Ixx,Iyy,Izz,de0(i),de1(i),de2(i),de3(i),e0(i),e1(i),e2(i),e3(i),g,h,m);
    adot_quat = ConsJacobdot_quat(de0(i),de1(i),de2(i),de3(i));
    a_quat = ConsJacob_quat(e0(i),e1(i),e2(i),e3(i));
    Q_quat = GenForc_quat(Moment,e0(i),e1(i),e2(i),e3(i));
    % We use Augmented Method to solve the dynamics of the system
    F = Q_quat - B_quat;
    A_Au = [M_quat -a_quat.'; -a_quat 0];
    B_Au = [F; adot_quat*dq_quat];
    dde = A_Au\B_Au;
    Endot_quat(i,1) = Energydot_quat(Ixx,Iyy,Izz,dde(1),dde(2),dde(3),dde(4),de0(i),de1(i),de2(i),de3(i),e0(i),e1(i),e2(i),e3(i),g,h,m);
end
P_quat = Power_quat(Moment,de0,de1,de2,de3,e0,e1,e2,e3);
figure;
plot(t,(Endot_quat-P_quat),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Power error','fontsize',14,'fontweight','bold');
end
end



%_____________________Modified Lagrange Method_____________________________
if ML_quat == 1
z0 = [1;0;0;0;0;-2;2];      %States: Z = [e0;e1;e2;e3;u1_quat;u2_quat;u3_quat]
options=odeset('maxstep',0.001);
[t,z] = ode45(@ModlagDyn_quat,simulation_time,z0, options);
e0_ML = z(:,1); e1_ML = z(:,2); e2_ML = z(:,3); e3_ML = z(:,4);
u1_quat = z(:,5); u2_quat = z(:,6); u3_quat = z(:,7);
de_Ml = zeros(4,length(simulation_time));
for i=1:length(simulation_time)
u_quat = [u1_quat(i);u2_quat(i);u3_quat(i)];
W_quat = W_Ml_quat(e0_ML(i),e1_ML(i),e2_ML(i),e3_ML(i));
X_quat = X_Ml_quat;
de_Ml(:,i) = W_quat*u_quat + X_quat;
end
%Constraint
f = zeros(length(simulation_time),1);
for i=1:length(simulation_time)
    a_quat = ConsJacob_quat(e0_ML(i),e1_ML(i),e2_ML(i),e3_ML(i));
    dq_quat = [de_Ml(1,i);de_Ml(2,i);de_Ml(3,i);de_Ml(4,i)];
    f(i,1) = a_quat*dq_quat;
end   
figure;
hold on;
plot(t,e0_ML ,'r-','linewidth',4)
plot(t,e1_ML ,'b-','linewidth',4)
plot(t,e2_ML ,'k-','linewidth',4)
plot(t,e3_ML ,'g-','linewidth',4)
legend('e_0','e_1','e_2','e_3');
set(gca,'fontsize',16,'fontweight','bold');
xlabel('Time(s)');
ylabel('Quaternions ','fontsize',14,'fontweight','bold');
figure;
plot(t,f,'k-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Error of Constraint','fontsize',14,'fontweight','bold');

if Moment == 0
En_Ml_quat = Energy_quat(Ixx,Iyy,Izz,de_Ml(1,:).',de_Ml(2,:).',de_Ml(3,:).',de_Ml(4,:).',e0_ML,e1_ML,e2_ML,e3_ML,g,h,m);
figure;
plot(t,100*(En_Ml_quat-En_Ml_quat(1))/En_Ml_quat(1),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Relative Error of  Mechanical Energy % ','fontsize',14,'fontweight','bold');
else
Endot_quat_itu = zeros(length(simulation_time),1);
for i = 1:length(simulation_time)
    u_quat = [u1_quat(i);u2_quat(i);u3_quat(i)];
    Qprime_ML_quat = Qpr_Ml_quat(Moment);
    Mprime_ML_quat = Mpr_Ml_quat(Ixx,Iyy,Izz,h,m); 
    Bprime_ML_quat_itu = Bpr_Ml_quat_itu(Ixx,Iyy,Izz,e0_ML(i),e1_ML(i),e2_ML(i),e3_ML(i),g,h,m,u1_quat(i),u2_quat(i),u3_quat(i));
    udot_quat = Mprime_ML_quat\(Qprime_ML_quat - Bprime_ML_quat_itu);
    Endot_quat_itu(i,1) = Energydot_quat_itu(Ixx,Iyy,Izz,udot_quat(1),udot_quat(2),udot_quat(3),e0_ML(i),e1_ML(i),e2_ML(i),e3_ML(i),g,h,m,u1_quat(i),u2_quat(i),u3_quat(i));
end
P_quat_itu = Power_quat_itu(Moment,u3_quat);
figure;
plot(t,(Endot_quat_itu-P_quat_itu),'r-','linewidth',1);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('Time (s)','fontsize',14,'fontweight','bold');
ylabel('Power error','fontsize',14,'fontweight','bold');
end
end
