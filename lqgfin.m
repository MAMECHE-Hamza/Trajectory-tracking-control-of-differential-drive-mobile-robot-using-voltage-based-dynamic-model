% MOBILE ROBOT PLANT SETUP AND LQG CONTROLLER DESIGN
clear all;
close all;
kx=0.1;
ky=0.1;
lx=0.3;
ly=0.1;
% Motors Specification
Kb1 = 0.71142; Kb2 = 0.70005;
Ra1 = 4; Ra2 = 3.5327;
Km1 = 0.95*Kb1; Km2 = 0.95*Kb2;
La1 = 0.0051381; La2 = 0.009607;
Bm1 = 0.028616; Bm2 = 0.044517;
Tf1 = 0.12256; Tf2 = 0.050745;
% Motor Driver lockup table
u_51=[[-5:0.2:-3],[-2.9:0.1:2.9],[3:0.2:5]];
vol_r=[-[10.31,9.98,9.88,9.8,9.72,9.62,9.5,9.37,9.22,9.05,8.87,8.73,...
    8.62,8.49,8.35,8.19,8.01,7.84,7.63,7.41,7.16,6.76,6.43,6.09,5.68,...
    5.22,4.74,4.16,3.56,2.92,2.15,0.42,0.35,0.032,0.027,0.023,0.019,...
    0.013,0.0068,0.00156],[0,0.0013,0.0018,0.006,0.012,0.018,0.026,...
    0.032,0.109,1.16,2,2.67,3.35,4.05,4.61,5.12,5.59,6.00,6.35,6.69,...
    7.06,7.29,7.55,7.77,7.96,8.14,8.31,8.45,8.59,8.72,8.86,9.09,9.25,...
    9.4,9.53,9.65,9.76,9.85,9.94,10.04,10.37]];
vol_l=[-[10.3,9.95,9.86,9.8,9.74,9.65,9.53,9.42,9.3,9.15,8.98,8.86,...
    8.75,8.64,8.53,8.39,8.25,8.05,7.9,7.6,7.4,7.1,6.88,6.48,6.1,5.7,...
    5.2,4.7,4.13,3.45,2.7,1.75,0.21,0.165,0.12,0.088,0.06,0.036,0.017,...
    0.005],[0,0.005,0.017,0.035,0.061,0.092,0.122,0.160,0.960,1.7,2.6,...
    3.25,4,4.6,5.1,5.6,6,6.44,6.8,7.1,7.4,7.65,7.8,8,8.11,8.25,8.4,...
    8.55,8.68,8.82,8.95,9.15,9.29,9.44,9.56,9.67,9.75,9.82,9.92,10.02,...
    10.36]];
% Robot Specification
r = 0.13/2;                %wheel radius in Meters
m = 8.2;                   %Mass in Kg
l=0.21;                    %Wheels axis length in Meters
I = 1.1114;                %Robot inertia
a=0.18;                    %Interest point distance from wheels axis 
% State-space model of the dynamic model 
A=[-Ra1/La1     0            -Kb1/(r*La1)          -Kb1*l/(2*r*La1)
    0           -Ra2/La2     -Kb2/(r*La2)           Kb2*l/(2*r*La2)
    Km1/(m*r)   Km2/(m*r)    -(Bm1+Bm2)/(m*r^2)   l*(Bm2-Bm1)/(2*m*r^2)
    Km1*l/(r*I) -Km2*l/(r*I) l*(Bm2-Bm1)/(I*r^2) -(l^2)*(Bm2+Bm1)/(2*I*r^2)];
B=[1/La1 0
   0     1/La2
   0     0
   0     0];
C=[0    0    1/r  l/(2*r)    
   0    0    1/r  -l/(2*r)];
D=[0 0
   0 0];
sys= ss(A,B,C,D);
sys.InputName={'RvoltageR','RvoltaeL'};
sys.OutputName={'OmegaMR','OmegaML'};
sys.StateName={'Ia1','Ia2','Lvelocity','Avelocity'};
% LQG parameters design
% LQR
WR_max=10;WL_max=10;       %rad/s
Ia1_max=1;Ia2_max=1;       %A
V_max=(WR_max+WL_max)*r/2;         %m/s
W_max=(WR_max+WL_max)*r/l;         %rad/s
alpha_ia1=0.1;alpha_ia2=0.1;alpha_v=0.1;alpha_w=0.5;
alpha_xi1=0.6;alpha_xi2=0.6;  
Vr_max=10;Vl_max=10;       %V
beta_vr=1/sqrt(2);beta_vl=1/sqrt(2);
rho=0.1;
v_q=[(alpha_ia1/Ia1_max)^2,(alpha_ia2/Ia2_max)^2,(alpha_v/V_max)^2,...
    (alpha_w/W_max)^2,(alpha_xi1)^2,(alpha_xi2)^2];
Q=diag(v_q);
v_r=[(beta_vr/Vr_max)^2,(beta_vl/Vl_max)^2];
R=rho*diag(v_r);
[K_lqi,S,e] = lqi(sys,Q,R);   
K_x=K_lqi(:,1:4);
% Kalman filter
sys_noise = ss(A,[B B],C,0);
sys_noise.InputName={'RvoltageR','RvoltaeL','volR_noise','volL_noise'};
sys_noise.OutputName={'OmegaMR','OmegaML'};
sys_noise.StateName={'Ia1','Ia2','Lvelocity','Avelocity'};
q1=0.1;q2=0.1;r1=0.1;r2=0.1;
Wn = [q1 0;0 q2];
Vn = [r1 0;0 r2];
[kalmf,L,P,M] = kalman(sys_noise,Wn,Vn);
[An,Bn,Cn,Dn] = ssdata(kalmf);
Cy=[1 0 0 0 0 0;0 1 0 0 0 0];      %extract y_e from kalmf output
Cx=[0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0; 0 0 0 0 0 1];    %extract x_e from kalmf output



       
   
     
    
    
