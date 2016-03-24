%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% This is solving a projectile problem by RK1, RK2, RK3, RK4 and RK5
close all; clear all; clc;
%% exact solution
syms x;
y0=200; x0=200; g=9.8; v0=20; gamma0=75*pi/180;
y(x)=y0+(x-x0)*tan(gamma0)-0.5*g/v0^2/(cos(gamma0))^2*(x-x0)^2;
K(x)=diff(diff(y(x)))/(1+(diff(y(x))).^2).^1.5;
v(x)=diff(y(x));
Gamma(x)=acos(v(x).^2/K(x)/g);
%% numerical solution parameter
A=[0,0,0,0;1,0,0,0;0,0,0,0;0,0,1,0];
B=[0;0;-g;0];
X0=[v0*cos(gamma0);x0;v0*sin(gamma0);y0];
order=4;
t_initial=0;
t_final=20;
dt=1e-1;
%% RK1
sol1(:,1)=X0;
m1=2;
while sol1(4,m1-1) >0
    [ S1, t1 ] = RK1( A,B,sol1(:,m1-1),dt,dt*(m1-1),dt*m1,order );
    sol1(:,m1)=S1(:,2);
    m1=m1+1;
end
figure(1);
set(gcf,'Color','w')
hold all;
plot(sol1(2,:),sol1(4,:),'red');
grid on;
for nn1=1:m1-1
    y_exact1(nn1)=double(y(sol1(2,nn1)));
end
plot(sol1(2,:),y_exact1,'blue')
legend('RK1 Solution','Exact Solution')
title('RK1 Solution','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('Y','Fontsize',20);
%% RK2
sol2(:,1)=X0;
m2=2;
while sol2(4,m2-1) >0
    [ S2, t2 ] = RK2( A,B,sol2(:,m2-1),dt,dt*(m2-1),dt*m2,order );
    sol2(:,m2)=S2(:,2);
    m2=m2+1;
end
figure(2);
set(gcf,'Color','w')
hold all;
plot(sol2(2,:),sol2(4,:),'red');
grid on;
for nn2=1:m2-1
    y_exact2(nn2)=double(y(sol2(2,nn2)));
end
plot(sol2(2,:),y_exact2,'blue')
legend('RK2 Solution','Exact Solution')
title('RK2 Solution','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('Y','Fontsize',20);
%% RK3
sol3(:,1)=X0;
m3=2;
while sol3(4,m3-1) >0
    [ S3, t3 ] = RK3( A,B,sol3(:,m3-1),dt,dt*(m3-1),dt*m3,order );
    sol3(:,m3)=S3(:,2);
    m3=m3+1;
end
figure(3);
set(gcf,'Color','w')
hold all;
plot(sol3(2,:),sol3(4,:),'red');
grid on;
for nn3=1:m3-1
    y_exact3(nn3)=double(y(sol3(2,nn3)));
end
plot(sol3(2,:),y_exact3,'blue')
legend('RK3 Solution','Exact Solution')
title('RK3 Solution','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('Y','Fontsize',20);
%% RK4
sol4(:,1)=X0;
m4=2;
while sol4(4,m4-1) >0
    [ S4, t4 ] = RK4( A,B,sol4(:,m4-1),dt,dt*(m4-1),dt*m4,order );
    sol4(:,m4)=S4(:,2);
    m4=m4+1;
end
figure(4);
set(gcf,'Color','w')
hold all;
plot(sol4(2,:),sol4(4,:),'red');
grid on;
for nn4=1:m4-1
    y_exact4(nn4)=double(y(sol4(2,nn4)));
end
plot(sol4(2,:),y_exact4,'blue')
legend('RK4 Solution','Exact Solution')
title('RK4 Solution','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('Y','Fontsize',20);
%% RK5
sol5(:,1)=X0;
m5=2;
while sol5(4,m5-1) >0
    [ S5, t5 ] = RK5( A,B,sol5(:,m5-1),dt,dt*(m5-1),dt*m5,order );
    sol5(:,m5)=S5(:,2);
    m5=m5+1;
end
figure(5);
set(gcf,'Color','w')
hold all;
plot(sol5(2,:),sol5(4,:),'red');
grid on;
for nn5=1:m5-1
    y_exact5(nn5)=double(y(sol5(2,nn5)));
end
plot(sol5(2,:),y_exact5,'blue')
legend('RK5 Solution','Exact Solution')
title('RK5 Solution','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('Y','Fontsize',20);
%% error
% RK1
[ e1, Max_e1, std_e1, mean_e1, RMS_e1, e_abs1, Max_abs_e1, std_abs_e1, mean_abs_e1, RMS_abs_e1 ] = ERROR ( sol1(4,:), y_exact1 );
figure(6);
set(gcf,'Color','w')
plot(sol1(2,:),e1);
legend(['Max. error = ' num2str(Max_e1) ', Std. of error = ' num2str(std_e1) ', Error mean = ' num2str(mean_e1) ', Error RMS = ' num2str(RMS_e1)])
title('RK1 Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('RK1_s_o_l - Exact_s_o_l','Fontsize',20);
grid on;
figure(7);
set(gcf,'Color','w')
plot(sol1(2,:),e_abs1);
legend(['Max. |error| = ' num2str(Max_abs_e1) ', Std. of |error| = ' num2str(std_abs_e1) ', Error mean = ' num2str(mean_abs_e1) ', Error RMS = ' num2str(RMS_abs_e1)])
title('RK1 absolute Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('|RK1_s_o_l - Exact_s_o_l|','Fontsize',20);
grid on;
% RK2
[ e2, Max_e2, std_e2, mean_e2, RMS_e2, e_abs2, Max_abs_e2, std_abs_e2, mean_abs_e2, RMS_abs_e2 ] = ERROR ( sol2(4,:), y_exact2 );
figure(8);
set(gcf,'Color','w')
plot(sol2(2,:),e2);
legend(['Max. error = ' num2str(Max_e2) ', Std. of error = ' num2str(std_e2) ', Error mean = ' num2str(mean_e2) ', Error RMS = ' num2str(RMS_e2)])
title('RK2 Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('RK2_s_o_l - Exact_s_o_l','Fontsize',20);
grid on;
figure(9);
set(gcf,'Color','w')
plot(sol2(2,:),e_abs2);
legend(['Max. |error| = ' num2str(Max_abs_e2) ', Std. of |error| = ' num2str(std_abs_e2) ', Error mean = ' num2str(mean_abs_e2) ', Error RMS = ' num2str(RMS_abs_e2)])
title('RK2 absolute Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('|RK2_s_o_l - Exact_s_o_l|','Fontsize',20);
grid on;
% RK3
[ e3, Max_e3, std_e3, mean_e3, RMS_e3, e_abs3, Max_abs_e3, std_abs_e3, mean_abs_e3, RMS_abs_e3 ] = ERROR ( sol3(4,:), y_exact3 );
figure(10);
set(gcf,'Color','w')
plot(sol3(2,:),e3);
legend(['Max. error = ' num2str(Max_e3) ', Std. of error = ' num2str(std_e3) ', Error mean = ' num2str(mean_e3) ', Error RMS = ' num2str(RMS_e3)])
title('RK3 Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('RK3_s_o_l - Exact_s_o_l','Fontsize',20);
grid on;
figure(11);
set(gcf,'Color','w')
plot(sol3(2,:),e_abs3);
legend(['Max. |error| = ' num2str(Max_abs_e3) ', Std. of |error| = ' num2str(std_abs_e3) ', Error mean = ' num2str(mean_abs_e3) ', Error RMS = ' num2str(RMS_abs_e3)])
title('RK3 absolute Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('|RK3_s_o_l - Exact_s_o_l|','Fontsize',20);
grid on;
% RK4
[ e4, Max_e4, std_e4, mean_e4, RMS_e4, e_abs4, Max_abs_e4, std_abs_e4, mean_abs_e4, RMS_abs_e4 ] = ERROR ( sol4(4,:), y_exact4 );
figure(12);
set(gcf,'Color','w')
plot(sol4(2,:),e4);
legend(['Max. error = ' num2str(Max_e4) ', Std. of error = ' num2str(std_e4) ', Error mean = ' num2str(mean_e4) ', Error RMS = ' num2str(RMS_e4)])
title('RK4 Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('RK4_s_o_l - Exact_s_o_l','Fontsize',20);
grid on;
figure(13);
set(gcf,'Color','w')
plot(sol4(2,:),e_abs4);
legend(['Max. |error| = ' num2str(Max_abs_e4) ', Std. of |error| = ' num2str(std_abs_e4) ', Error mean = ' num2str(mean_abs_e4) ', Error RMS = ' num2str(RMS_abs_e4)])
title('RK4 absolute Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('|RK4_s_o_l - Exact_s_o_l|','Fontsize',20);
grid on;
% RK5
[ e5, Max_e5, std_e5, mean_e5, RMS_e5, e_abs5, Max_abs_e5, std_abs_e5, mean_abs_e5, RMS_abs_e5 ] = ERROR ( sol5(4,:), y_exact5 );
figure(14);
set(gcf,'Color','w')
plot(sol5(2,:),e5);
legend(['Max. error = ' num2str(Max_e5) ', Std. of error = ' num2str(std_e5) ', Error mean = ' num2str(mean_e5) ', Error RMS = ' num2str(RMS_e5)])
title('RK5 Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('RK5_s_o_l - Exact_s_o_l','Fontsize',20);
grid on;
figure(15);
set(gcf,'Color','w')
plot(sol5(2,:),e_abs5);
legend(['Max. |error| = ' num2str(Max_abs_e5) ', Std. of |error| = ' num2str(std_abs_e5) ', Error mean = ' num2str(mean_abs_e5) ', Error RMS = ' num2str(RMS_abs_e5)])
title('RK5 absolute Error','Fontsize',18);
xlabel('X','Fontsize',20);
ylabel('|RK5_s_o_l - Exact_s_o_l|','Fontsize',20);
grid on;











