syms theta h eta phi r alpha beta gamma
%% variables

% 1 

phi = pi/10;
eta = pi/10;
theta= pi/10;

% 2

alpha = pi/20;
beta = pi/20;
gamma = pi/20;

% 2

alpha3 = 0;
beta3 = 0;
gamma3 = 0;

% lengths

r = 1;
h = 1;

%% ## BASE COORDINATES ##

Ori = [0 0 0];

b0 = [0;0;-h]*r;
b1 = [cos(0*pi/2);sin(0*pi/2);-h]*r;
b2 = [cos(1*pi/2);sin(1*pi/2);-h]*r;
b3 = [cos(2*pi/2);sin(2*pi/2);-h]*r;
b4 = [cos(3*pi/2);sin(3*pi/2);-h]*r;

Rotz =[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1] ;

B0 = Rotz*b0;
B1 = Rotz*b1;
B2 = Rotz*b2;
B3 = Rotz*b3;
B4 = Rotz*b4 ;

%% Plane 2 COORDINATES ##

a1 = [cos(0*pi/2);sin(0*pi/2);h]*r;
a2 = [cos(1*pi/2);sin(1*pi/2);h]*r;
a3 = [cos(2*pi/2);sin(2*pi/2);h]*r;
a4 = [cos(3*pi/2);sin(3*pi/2);h]*r;
a0 = [0;0;h]*r;

A1 = a1;
A2 = a2;
A3 = a3;
A4 = a4;
A0 = a0;

rotX = [1, 0,0;0, cos(eta), -sin(eta);0,sin(eta), cos(eta)] ;
rotY = [cos(phi), 0, sin(phi);0,1,0;-sin(phi), 0,cos(phi)] ;
rotXY = rotX*rotY ;

C1 =  Rotz*rotXY*A1 ;
C2 =  Rotz*rotXY*A2 ;
C3 =  Rotz*rotXY*A3 ;
C4 =  Rotz*rotXY*A4 ;
C0 =  Rotz*rotXY*A0 ;

%% Plane 3 COORDINATES %%
P1 =  [cos(0*pi/2);sin(0*pi/2);3*h]*r ;
P2 =  [cos(1*pi/2);sin(1*pi/2);3*h]*r ;
P3 =  [cos(2*pi/2);sin(2*pi/2);3*h]*r ;
P4 =  [cos(3*pi/2);sin(3*pi/2);3*h]*r ;
P0 =  [0;0;3*h]*r ;
Ori_2 = [0;0;2*h]*r;


Rotz2 =[cos(gamma),-sin(gamma),0;sin(gamma),cos(gamma),0;0,0,1] ;
rotX2 = [1, 0,0;0, cos(alpha), -sin(alpha);0,sin(alpha), cos(alpha)] ;
rotY2 = [cos(beta), 0, sin(beta);0,1,0;-sin(beta), 0,cos(beta)] ;
rotXY2 = rotX2*rotY2 ;

Q1 =  Rotz2*rotXY2*P1 ;
Q2 =  Rotz2*rotXY2*P2 ;
Q3 =  Rotz2*rotXY2*P3 ;
Q4 =  Rotz2*rotXY2*P4 ;
Q0 =  Rotz2*rotXY2*P0 ;

Ori2 = Rotz2*rotXY2*Ori_2 ;

%% Plane 4 COORDINATES %%
R1 =  [cos(0*pi/2);sin(0*pi/2);5*h]*r ;
R2 =  [cos(1*pi/2);sin(1*pi/2);5*h]*r ;
R3 =  [cos(2*pi/2);sin(2*pi/2);5*h]*r ;
R4 =  [cos(3*pi/2);sin(3*pi/2);5*h]*r ;
R0 =  [0;0;5*h]*r ;
Ori_3 = [0;0;4*h]*r;


Rotz3 =[cos(gamma3),-sin(gamma3),0;sin(gamma3),cos(gamma3),0;0,0,1] ;
rotX3 = [1, 0,0;0, cos(alpha3), -sin(alpha3);0,sin(alpha3), cos(alpha3)] ;
rotY3 = [cos(beta3), 0, sin(beta3);0,1,0;-sin(beta3), 0,cos(beta3)] ;
rotXY3 = rotX3*rotY3 ;

S1 =  Rotz3*rotXY3*R1 ;
S2 =  Rotz3*rotXY3*R2 ;
S3 =  Rotz3*rotXY3*R3 ;
S4 =  Rotz3*rotXY3*R4 ;
S0 =  Rotz3*rotXY3*R0 ;

Ori3 = Rotz3*rotXY3*Ori_3 ;


%% IKP %%

% 1

L1 =  (sqrt((B1(1)-C1(1))^2+(B1(2)-C1(2))^2+(B1(3)-C1(3))^2)) ;
L2 =  (sqrt((B2(1)-C2(1))^2+(B2(2)-C2(2))^2+(B2(3)-C2(3))^2)) ;
L3 =  (sqrt((B3(1)-C3(1))^2+(B3(2)-C3(2))^2+(B3(3)-C3(3))^2)) ;
L4 =  (sqrt((B4(1)-C4(1))^2+(B4(2)-C4(2))^2+(B4(3)-C4(3))^2)) ;

% 2

L5 =  sqrt((C1(1)-Q1(1))^2+(C1(2)-Q1(2))^2+(C1(3)-Q1(3))^2) ;
L6 =  sqrt((C2(1)-Q2(1))^2+(C2(2)-Q2(2))^2+(C2(3)-Q2(3))^2) ;
L7 =  sqrt((C3(1)-Q3(1))^2+(C3(2)-Q3(2))^2+(C3(3)-Q3(3))^2) ;
L8 =  sqrt((C4(1)-Q4(1))^2+(C4(2)-Q4(2))^2+(C4(3)-Q4(3))^2) ;

% 2

L9 =   sqrt((Q1(1)-S1(1))^2+(Q1(2)-S1(2))^2+(Q1(3)-S1(3))^2) ;
L10 =  sqrt((Q2(1)-S2(1))^2+(Q2(2)-S2(2))^2+(Q2(3)-S2(3))^2) ;
L11 =  sqrt((Q3(1)-S3(1))^2+(Q3(2)-S3(2))^2+(Q3(3)-S3(3))^2) ;
L12 =  sqrt((Q4(1)-S4(1))^2+(Q4(2)-S4(2))^2+(Q4(3)-S4(3))^2) ;

 %% plotting
light('style','local','position',[0 0 3],'color',[1 1 1])
axis vis3d
box on 
grid off
%% plot joints
facecolor='none';
    edgecolor='k';
  jointfc='k';jointfc1='b';
   jointec='k';
linkfc='r';
 linkfc2='r';
linkec='k';
linkec1='none';

k1_A1=B1'-Ori/norm(B1'-Ori);
k2_A2=B2'-Ori/norm(B2'-Ori);
k3_A3=B3'-Ori/norm(B3'-Ori);
k4_A4=B4'-Ori/norm(B4'-Ori);
k0_A0=B0'-Ori/norm(B0'-Ori);


k1_C1=C1'-Ori/norm(C1'-Ori);
k2_C2=C2'-Ori/norm(C2'-Ori);
k3_C3=C3'-Ori/norm(C3'-Ori);
k4_C4=C4'-Ori/norm(C4'-Ori);
k0_C0=C0'-Ori/norm(C0'-Ori);

k1_Q1=Q1'-Ori/norm(Q1'-Ori);
k2_Q2=Q2'-Ori/norm(Q2'-Ori);
k3_Q3=Q3'-Ori/norm(Q3'-Ori);
k4_Q4=Q4'-Ori/norm(Q4'-Ori);
k0_Q0=Q0'-Ori/norm(Q0'-Ori);
k0_Ori2=Ori2'-Ori/norm(Ori2'-Ori);

k1_S1=S1'-Ori/norm(S1'-Ori);
k2_S2=S2'-Ori/norm(S2'-Ori);
k3_S3=S3'-Ori/norm(S3'-Ori);
k4_S4=S4'-Ori/norm(S4'-Ori);
k0_S0=S0'-Ori/norm(S0'-Ori);
k0_Ori3=Ori3'-Ori/norm(Ori3'-Ori);

%% joints on the ground plane

sphere_joint_axis(10,k1_A1,0.03,eye(3),Ori,jointfc,jointec)
hold on

sphere_joint_axis(10,k1_A1,0.03,eye(3),B1',jointfc1,jointec)
hold on

sphere_joint_axis(10,k2_A2,0.03,eye(3),B2',jointfc1,jointec)
hold on

sphere_joint_axis(10,k3_A3,0.03,eye(3),B3',jointfc1,jointec)
hold on

sphere_joint_axis(10,k4_A4,0.03,eye(3),B4',jointfc1,jointec)
hold on

sphere_joint_axis(10,k0_A0,0.03,eye(3),B0',jointfc1,jointec)
hold on

%% joints on the 1 plane

sphere_joint_axis(10,k1_C1,0.03,eye(3),C1',jointfc1,jointec)
hold on

sphere_joint_axis(10,k2_C2,0.03,eye(3),C2',jointfc1,jointec)
hold on

sphere_joint_axis(10,k3_C3,0.03,eye(3),C3',jointfc1,jointec)
hold on

sphere_joint_axis(10,k4_C4,0.03,eye(3),C4',jointfc1,jointec)
hold on

sphere_joint_axis(10,k0_C0,0.03,eye(3),C0',jointfc1,jointec)
hold on


%% joints on the 2 plane

sphere_joint_axis(10,k1_Q1,0.03,eye(3),Q1',jointfc1,jointec)
hold on

sphere_joint_axis(10,k2_Q2,0.03,eye(3),Q2',jointfc1,jointec)
hold on

sphere_joint_axis(10,k3_Q3,0.03,eye(3),Q3',jointfc1,jointec)
hold on

sphere_joint_axis(10,k4_Q4,0.03,eye(3),Q4',jointfc1,jointec)
hold on

sphere_joint_axis(10,k0_Q0,0.03,eye(3),Q0',jointfc1,jointec)
hold on

sphere_joint_axis(10,k0_Ori2,0.03,eye(3),Ori2',jointfc1,jointec)
hold on

%% joints on the 3 plane

sphere_joint_axis(10,k1_S1,0.03,eye(3),S1',jointfc1,jointec)
hold on

sphere_joint_axis(10,k2_S2,0.03,eye(3),S2',jointfc1,jointec)
hold on

sphere_joint_axis(10,k3_S3,0.03,eye(3),S3',jointfc1,jointec)
hold on

sphere_joint_axis(10,k4_S4,0.03,eye(3),S4',jointfc1,jointec)
hold on

sphere_joint_axis(10,k0_S0,0.03,eye(3),S0',jointfc1,jointec)
hold on

sphere_joint_axis(10,k0_Ori3,0.03,eye(3),Ori3',jointfc1,jointec)
hold on
%% connections

cylinderbetweenpoints(0.03,10,[B0(1),B0(2),B0(3)],[Ori(1),Ori(2),Ori(3)],'r',linkec1);
hold on
cylinderbetweenpoints(0.03,10,[C0(1),C0(2),C0(3)],[Ori(1),Ori(2),Ori(3)],'r',linkec1);
hold on

cylinderbetweenpoints(0.03,10,[Ori2(1),Ori2(2),Ori2(3)],[Q0(1),Q0(2),Q0(3)],'r',linkec1);
hold on
cylinderbetweenpoints(0.03,10,[C0(1),C0(2),C0(3)],[Ori2(1),Ori2(2),Ori2(3)],'r',linkec1);
hold on


cylinderbetweenpoints(0.03,10,[Ori3(1),Ori3(2),Ori3(3)],[S0(1),S0(2),S0(3)],'r',linkec1);
hold on
cylinderbetweenpoints(0.03,10,[Q0(1),Q0(2),Q0(3)],[Ori3(1),Ori3(2),Ori3(3)],'r',linkec1);
hold on

%% notations

text(B1(1)+0.03,B1(2)+0.03,B1(3),'B1','fontsize',14,'Interpreter', 'latex');
text(B2(1)+0.03,B2(2)+0.03,B2(3),'B2','fontsize',14,'Interpreter', 'latex');
text(B3(1)+0.03,B3(2)+0.03,B3(3),'B3','fontsize',14,'Interpreter', 'latex');
text(B4(1)+0.03,B4(2)+0.03,B4(3),'B4','fontsize',14,'Interpreter', 'latex');
text(B0(1)+0.03,B0(2)+0.03,B0(3),'B0','fontsize',14,'Interpreter', 'latex');

text(C1(1)+0.03,C1(2)+0.03,C1(3),'C1','fontsize',14,'Interpreter', 'latex');
text(C2(1)+0.03,C2(2)+0.03,C2(3),'C2','fontsize',14,'Interpreter', 'latex');
text(C3(1)+0.03,C3(2)+0.03,C3(3),'C3','fontsize',14,'Interpreter', 'latex');
text(C4(1)+0.03,C4(2)+0.03,C4(3),'C4','fontsize',14,'Interpreter', 'latex');
text(C0(1)+0.03,C0(2)+0.03,C0(3),'C0','fontsize',14,'Interpreter', 'latex');

text(Q1(1)+0.03,Q1(2)+0.03,Q1(3),'Q1','fontsize',14,'Interpreter', 'latex');
text(Q2(1)+0.03,Q2(2)+0.03,Q2(3),'Q2','fontsize',14,'Interpreter', 'latex');
text(Q3(1)+0.03,Q3(2)+0.03,Q3(3),'Q3','fontsize',14,'Interpreter', 'latex');
text(Q4(1)+0.03,Q4(2)+0.03,Q4(3),'Q4','fontsize',14,'Interpreter', 'latex');
text(Q0(1)+0.03,Q0(2)+0.03,Q0(3),'Q0','fontsize',14,'Interpreter', 'latex');

text(S1(1)+0.03,S1(2)+0.03,S1(3),'S1','fontsize',14,'Interpreter', 'latex');
text(S2(1)+0.03,S2(2)+0.03,S2(3),'S2','fontsize',14,'Interpreter', 'latex');
text(S3(1)+0.03,S3(2)+0.03,S3(3),'S3','fontsize',14,'Interpreter', 'latex');
text(S4(1)+0.03,S4(2)+0.03,S4(3),'S4','fontsize',14,'Interpreter', 'latex');
text(S0(1)+0.03,S0(2)+0.03,S0(3),'S0','fontsize',14,'Interpreter', 'latex');
%% prismatic joint
hold on
prismatic_joint_axis(0.01,[B1(1),B1(2),B1(3)],[C1(1),C1(2),C1(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[B2(1),B2(2),B2(3)],[C2(1),C2(2),C2(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[B3(1),B3(2),B3(3)],[C3(1),C3(2),C3(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[B4(1),B4(2),B4(3)],[C4(1),C4(2),C4(3)],jointfc,jointec,linkfc)

%%
hold on
prismatic_joint_axis(0.01,[C1(1),C1(2),C1(3)],[Q1(1),Q1(2),Q1(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[C2(1),C2(2),C2(3)],[Q2(1),Q2(2),Q2(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[C3(1),C3(2),C3(3)],[Q3(1),Q3(2),Q3(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[C4(1),C4(2),C4(3)],[Q4(1),Q4(2),Q4(3)],jointfc,jointec,linkfc)

%%
hold on
prismatic_joint_axis(0.01,[Q1(1),Q1(2),Q1(3)],[S1(1),S1(2),S1(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[Q2(1),Q2(2),Q2(3)],[S2(1),S2(2),S2(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[Q3(1),Q3(2),Q3(3)],[S3(1),S3(2),S3(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[Q4(1),Q4(2),Q4(3)],[S4(1),S4(2),S4(3)],jointfc,jointec,linkfc)

%% surface

shadow=0;%1=shadow, 0=none
xa= [B1(1) B2(1) B3(1) B4(1)];
ya= [B1(2) B2(2) B3(2) B4(2)];
za= [B1(3) B2(3) B3(3) B4(3)];
SS= fill3(xa,ya,za,'r');
set(SS,'facealpha',1)


shadow=0;%1=shadow, 0=none
xa= [C1(1) C2(1) C3(1) C4(1)];
ya= [C1(2) C2(2) C3(2) C4(2)];
za= [C1(3) C2(3) C3(3) C4(3)];
SS= fill3(xa,ya,za,'r');
set(SS,'facealpha',1)

shadow=0;%1=shadow, 0=none
xa= [Q1(1) Q2(1) Q3(1) Q4(1)];
ya= [Q1(2) Q2(2) Q3(2) Q4(2)];
za= [Q1(3) Q2(3) Q3(3) Q4(3)];
SS= fill3(xa,ya,za,'r');
set(SS,'facealpha',1)

shadow=0;%1=shadow, 0=none
xa= [S1(1) S2(1) S3(1) S4(1)];
ya= [S1(2) S2(2) S3(2) S4(2)];
za= [S1(3) S2(3) S3(3) S4(3)];
SS= fill3(xa,ya,za,'r');
set(SS,'facealpha',1)


%%

