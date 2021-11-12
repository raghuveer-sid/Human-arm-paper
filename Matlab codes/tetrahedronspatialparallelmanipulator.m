    %% based on 2spsUrizar 2014 paper
% syms beta2 beta3 gamma3 phi theta psi real
% set(gcf, 'color', 'white','units','pixels','position',[0 0 1920 1080]);
%% first case
beta2=deg2rad(30);
beta3=deg2rad(30);
gamma3=deg2rad(60);
R=1;
L=0.5;
% beta2=deg2rad(0);
% gamma3=deg2rad(30);
% beta3=deg2rad(0);

% vidObj=VideoWriter('tetrahedronmanipulator.mp4');
% vidObj.FrameRate = 24;
% open(vidObj);
% %% general case
% beta2=deg2rad(30);
% beta3=deg2rad(30);
% gamma3=deg2rad(60);

%  for psi=-pi/6:0.01:pi/6
%   theta=deg2rad(60);
%  phi=deg2rad(30);
%  psi=deg2rad(30);
 theta=deg2rad(0);
 phi=deg2rad(0);
 psi=deg2rad(0);
% R=1;
% L=1;
OA_1=R;
OA_2=R;
OA_3=R;
Ori=[0,0,0];
 a1=[0,0,R]';
% a1=[-R,0,0]';
a2=[R,0,0]';
a3=[0,R,0]';
%% wrt moving frame uvw
b2u=cos(beta2);
b2v=0;
b2w=sin(beta2);
b3u=cos(beta3)*cos(gamma3);
b3v=cos(beta3)*sin(gamma3);
b3w=sin(beta3);
Mb1=L*[0,0,1]';
Mb2=L*[b2u,b2v,b2w]';
Mb3=L*[b3u,b3v,b3w]';

%% wvw euler angle zyz format (phi, theta, psi)
rotz_phi= [cos(phi) -sin(phi) 0
    sin(phi) cos(phi) 0
    0 0 1];
roty_theta= [cos(theta) 0 sin(theta)
    0 1 0
    -sin(theta) 0 cos(theta)];
rotz_psi= [cos(psi) -sin(psi) 0
    sin(psi) cos(psi) 0
    0 0 1];
R_zyz=rotz_phi*roty_theta*rotz_psi;
%% position vector of B expressed in fixed frame A
b1=R_zyz*Mb1;
b2=R_zyz*Mb2;
b3=R_zyz*Mb3;

%% inverse kinematics loop closure equation   euler angles theta phi psi known find limb lengths l1,l2,l3
% l1=sqrt(a1.^2+b1.^2-2*a1.'*b1);
% l2=sqrt(a2.^2+b2.^2-2*a2.'*b2);
% l3=sqrt(a3.^2+b3.^2-2*a3.'*b3);

l1=sqrt((b1-a1)'*(b1-a1));
l2=sqrt((b2-a2)'*(b2-a2));
l3=sqrt((b3-a3)'*(b3-a3));
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
k1_A1=a1'-Ori/norm(a1'-Ori)
k2_A2=a2'-Ori/norm(a2'-Ori);
k3_A3=a3'-Ori/norm(a3'-Ori);
k1_B1=b1'-Ori/norm(b1'-Ori);
k2_B2=b2'-Ori/norm(b2'-Ori);
k3_B3=b3'-Ori/norm(b3'-Ori);
%% joints on the ground plane
sphere_joint_axis(10,k1_A1,0.03,eye(3),Ori,jointfc,jointec)
hold on

sphere_joint_axis(10,k1_A1,0.03,eye(3),a1',jointfc1,jointec)
hold on

sphere_joint_axis(10,k2_A2,0.03,eye(3),a2',jointfc1,jointec)
hold on

sphere_joint_axis(10,k3_A3,0.03,eye(3),a3',jointfc1,jointec)
hold on

sphere_joint_axis(10,k1_B1,0.03,eye(3),b1',jointfc1,jointec)
hold on

sphere_joint_axis(10,k2_B2,0.03,eye(3),b2',jointfc1,jointec)
hold on

sphere_joint_axis(10,k3_B3,0.03,eye(3),b3',jointfc1,jointec)
hold on

cylinderbetweenpoints(0.03,10,[b1(1),b1(2),b1(3)],[b2(1),b2(2),b2(3)],'r',linkec1);
hold on
cylinderbetweenpoints(0.03,10,[b3(1),b3(2),b3(3)],[b2(1),b2(2),b2(3)],'r',linkec1);
hold on
cylinderbetweenpoints(0.03,10,[b1(1),b1(2),b1(3)],[b3(1),b3(2),b3(3)],'r',linkec1);
hold on
cylinderbetweenpoints(0.03,10,[b1(1),b1(2),b1(3)],[Ori(1),Ori(2),Ori(3)],'r',linkec);
hold on
cylinderbetweenpoints(0.03,10,[b2(1),b2(2),b2(3)],[Ori(1),Ori(2),Ori(3)],'r',linkec);
hold on
cylinderbetweenpoints(0.03,10,[b3(1),b3(2),b3(3)],[Ori(1),Ori(2),Ori(3)],'r',linkec);

text(b1(1)+0.03,b1(2)+0.03,b1(3),'b1','fontsize',14,'Interpreter', 'latex');
text(b2(1)+0.03,b2(2)+0.03,b2(3),'b2','fontsize',14,'Interpreter', 'latex');
text(b3(1)+0.03,b3(2)+0.03,b3(3),'b3','fontsize',14,'Interpreter', 'latex');

text(a1(1)+0.03,a1(2)+0.03,a1(3),'a1','fontsize',14,'Interpreter', 'latex');
text(a2(1)+0.03,a2(2)+0.03,a2(3),'a2','fontsize',14,'Interpreter', 'latex');
text(a3(1)+0.03,a3(2)+0.03,a3(3),'a3','fontsize',14,'Interpreter', 'latex');
%% prismatic joint
hold on
prismatic_joint_axis(0.01,[a1(1),a1(2),a1(3)],[b1(1),b1(2),b1(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[a2(1),a2(2),a2(3)],[b2(1),b2(2),b2(3)],jointfc,jointec,linkfc)

hold on
prismatic_joint_axis(0.01,[a3(1),a3(2),a3(3)],[b3(1),b3(2),b3(3)],jointfc,jointec,linkfc)


shadow=0;%1=shadow, 0=none
xa= [b1(1) b2(1) b3(1)];
ya= [b1(2) b2(2) b3(2)];
za= [b1(3) b2(3) b3(3)];
h= fill3(xa,ya,za,'r');
set(h,'facealpha',.5)
% xlabel('i')
% ylabel('j')
% zlabel('k')
set(gca, 'FontName', 'Times', 'FontSize', 18)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$z$', 'Interpreter', 'latex')
%Camera normal vector
    V1=xa;
    V2=ya;
    normal=cross(V1,V2)/norm(cross(V1,V2));
%     set(gca,'CameraUpVector',normal,'Projection','perspective','CameraViewAngle',30)
%             camoffset=4*(V1+V2);
%             campos(normal/1.5-camoffset)
%     %Setting camera target
%             target=(xa+ya)/2-normal;
%             camtarget([target])
%             camzoom(1)
%    view([-1.029843538673483e+02,21.415023228855262])  
%   view([1.560065460485506e+02,25.230548757282591])
%             view([1,1,1])
% set(gcf, 'Position',  [100, 100, 500, 400])
%             set(gca,'Projection','perspective','CameraViewAngle',30)
%             position=9*[.1 -1.3 .6];
%             campos(position)
%             camtarget([0 0 .7])
            camzoom(1.5)
%Plotting shadow
switch shadow
  case 1
shadowcolor=.95*[1 1 1];
patch(xa,ya,shadowcolor,'edgecolor','none')
otherwise
end
view([95.398134742922295,19.453786480907187])
set(gcf,'color','white')
drawnow
hold off
%   currFrame= getframe(gcf);
%           writeVideo(vidObj, currFrame);
%     hold on       
%  end
%       close(gcf)
%  close(vidObj);


%% direct kinematics
% given input variables and lengths, find output angles  the three euler
% angles
 m1= R^2+L^2-l1^2;
 theta=acos(m1/2*R*L);
 
 %% DKP Jaco
%  syms beta2 beta3 gamma3 phi theta psi R L real
% R=1; L=0.5;
R=1; L=1;
beta2=deg2rad(15);
gamma3=deg2rad(30);
beta3=deg2rad(15);
rotz_phi= [cos(phi) -sin(phi) 0
    sin(phi) cos(phi) 0
    0 0 1];
roty_theta= [cos(theta) 0 sin(theta)
    0 1 0
    -sin(theta) 0 cos(theta)];
rotz_psi= [cos(psi) -sin(psi) 0
    sin(psi) cos(psi) 0
    0 0 1];
R_zyz=rotz_phi*roty_theta*rotz_psi;
b2u=cos(beta2);
b2v=0;
b2w=sin(beta2);
b3u=cos(beta3)*cos(gamma3);
b3v=cos(beta3)*sin(gamma3);
b3w=sin(beta3);
Mb1=L*[0,0,1]';
Mb2=L*[b2u,b2v,b2w]';
Mb3=L*[b3u,b3v,b3w]';
b1=R_zyz*Mb1;
b2=R_zyz*Mb2;
b3=R_zyz*Mb3;
 a1=[0,0,R]';
a2=[R,0,0]';
a3=[0,R,0]';
 Jdkp=[cross(b1,(b1-a1))'
     cross(b2,(b2-a2))'
     cross(b3,(b3-a3))'];
 det(Jdkp)
Lde=simplify(subs(det(Jdkp)))
  simplify(det(Jdkp))
beta2=deg2rad(0);
beta3=deg2rad(0);
gamma3=deg2rad(90);

%% cusp case 
% R=1; L=0.5;
% beta2=deg2rad(10);
% gamma3=deg2rad(30);
% beta3=deg2rad(10);
% Lde=simplify(subs(det(Jdkp)));

figure(2)
Q=[0;0;0];

step=0.1;
% [phi,theta,psi] = meshgrid(linspace(-pi,step, pi), ...
%                          linspace(0,step,pi), ...
%                          linspace(-pi,step, pi));
for phi=-pi:step:pi
    for theta=pi/6
%     for theta=0:step:pi
        for psi=-pi:step:pi
            rotz_phi= [cos(phi) -sin(phi) 0;
    sin(phi) cos(phi) 0;
    0 0 1];
roty_theta= [cos(theta) 0 sin(theta);
    0 1 0;
    -sin(theta) 0 cos(theta)];
rotz_psi= [cos(psi) -sin(psi) 0;
    sin(psi) cos(psi) 0;
    0 0 1];
R_zyz=rotz_phi*roty_theta*rotz_psi;
b1=R_zyz*Mb1;
b2=R_zyz*Mb2;
b3=R_zyz*Mb3;
%             f2 =(cos(pi/18).*(cos(pi/18).*cos(phi).^2 - cos(pi/18).*cos(psi).^2 - cos(pi/18).*cos(phi).^2.*cos(theta).^2 + cos(pi/18).*cos(psi).^2.*cos(theta).^2 + sin(pi/18).*cos(psi).*cos(theta).*sin(theta) - 3^(1/2).*cos(pi/18).*cos(psi).*sin(psi) + sin(pi/18).*cos(phi).*sin(phi).*sin(psi).*sin(theta) - sin(pi/18).*cos(phi).^2.*cos(psi).*cos(theta).*sin(theta) + 2*sin(pi/18).*cos(phi).^2.*cos(theta).*sin(psi).*sin(theta) + 3^(1/2).*sin(pi/18).*cos(theta).*sin(psi).*sin(theta) + 3^(1/2)*cos(pi/18).*cos(psi).*cos(theta).^2.*sin(psi) + 2*sin(pi/18).*cos(phi).*cos(psi).*sin(phi).*sin(theta) - 3^(1/2).*sin(pi/18).*cos(phi).*cos(psi).*sin(phi).*sin(theta) - 3^(1/2)*sin(pi/18).*cos(phi).^2.*cos(theta).*sin(psi).*sin(theta)))/16;
            f2=-sin(theta).^2.*(sin(phi).^2 - sin(psi).^2);
if f2 ~= 0
%     Q1=[phi;theta;psi];
%     Q=[Q,Q1];
    l1=sqrt(((b1)-a1)'*((b1)-a1));
    l2=sqrt(((b2)-a2)'*((b2)-a2));
    l3=sqrt(((b3)-a3)'*((b3)-a3));
     Q1=[l1;l2;l3];
%      f=f2(phi,theta,psi);
     Q=[Q,Q1];
     
end
    
    
            %  Q=[0;0];
% for theta_a= -pi:step:pi
%     for theta_b= -pi:step:pi
%  if wrench_closure_combinatoric_null_space(double(subs(W)))
%      Q1=[theta_a;theta_b];
%      Q=[Q,Q1];
%  else
%      disp('wrench closure not satisfied')
%  end
%  
%     end
% end
% 
% c=length(Q);
% for i=1:1:c
%     Qx(i)=Q(1,i);
%     Qy(i)=Q(2,i);
%  
% end
% plot(Qx,Qy,'.r')
% xlim([-pi pi]); ylim([-pi pi]);



% [phi,theta,psi] = meshgrid(linspace(-pi, pi), ...
%                         linspace(0,pi), ...
%                         linspace(-pi, pi));
% f2 =@(phi,theta,psi)(cos(pi/18).*(cos(pi/18).*cos(phi).^2 - cos(pi/18).*cos(psi).^2 - cos(pi/18).*cos(phi).^2.*cos(theta).^2 + cos(pi/18).*cos(psi).^2.*cos(theta).^2 + sin(pi/18).*cos(psi).*cos(theta).*sin(theta) - 3^(1/2).*cos(pi/18).*cos(psi).*sin(psi) + sin(pi/18).*cos(phi).*sin(phi).*sin(psi).*sin(theta) - sin(pi/18).*cos(phi).^2.*cos(psi).*cos(theta).*sin(theta) + 2*sin(pi/18).*cos(phi).^2.*cos(theta).*sin(psi).*sin(theta) + 3^(1/2).*sin(pi/18).*cos(theta).*sin(psi).*sin(theta) + 3^(1/2)*cos(pi/18).*cos(psi).*cos(theta).^2.*sin(psi) + 2*sin(pi/18).*cos(phi).*cos(psi).*sin(phi).*sin(theta) - 3^(1/2).*sin(pi/18).*cos(phi).*cos(psi).*sin(phi).*sin(theta) - 3^(1/2)*sin(pi/18).*cos(phi).^2.*cos(theta).*sin(psi).*sin(theta)))/16;
% if ~isempty(f2)
%     l1=sqrt((b1-a1)'*(b1-a1));
%     l2=sqrt((b2-a2)'*(b2-a2));
%     l3=sqrt((b3-a3)'*(b3-a3));
% end
        end
    end
end
c=length(Q);
for i=1:1:c
    Qx(i)=Q(1,i);
    Qy(i)=Q(2,i);
    Qz(i)=Q(3,i);
end
plot3(Qx,Qy,Qz,'.r')
scatter3(Qx,Qy,Qz,'.r')
xlim([abs(R-L) R+L]); ylim([abs(R-L) R+L]);  zlim([abs(R-L) R+L]);
tri = delaunay(Qx,Qy);
% plot(Qx,Qy,'.')

%%
% How many triangles are there?

[r,c] = size(tri);
disp(r)

%% Plot it with TRISURF
C = 0.7.*ones(1,size(Qz,2));
 fill3(Qx,Qy,Qz,C) 
%  h = trisurf(tri, Qx,Qy, Qz);
% h = trisurf(tri, Qx,Qy, Qz,'Facecolor','red','FaceAlpha',.1,'EdgeColor','none');

% h.EdgeColor = 'none';
% h.FaceColor = 'interp';
% h.FaceLighting = 'gouraud';
% h.SpecularStrength = 5/8;
% light('Position',[9 -5 8])

axis vis3d
xlabel('l1(mm)')
ylabel('l2(mm)')
zlabel('l3(mm)')
% lighting phong
% shading interp

% hel3 = isosurface(Qx,Qy,Qz, f, 0);
% patch(hel3,'FaceColor',[1 .5 0],'EdgeColor','none');
% view(3)
% camlight
% ax = gca;
% set(ax,'XLim',[abs(R-L) R+L],'YLim',[abs(R-L) R+L],'ZLim',[abs(R-L) R+L],'DataAspectRatio',[1 1 1])
% box on
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('z(mm)')
% interval = [-pi pi 0 pi -pi pi];
%  fimplicit3(f2,interval,'EdgeColor','none','FaceAlpha',0.8)
% % fsurf(f2,interval,'EdgeColor','none','FaceAlpha',1)
% shading interp;
%     light;
%     lighting phong;
% hel1 = isosurface(phi,theta,psi,f2, 0);
% patch(hel1,'FaceColor',[.5 1 .5],'EdgeColor','none');
% 
% view(3)
% camlight
% ax = gca;
% set(ax,'XLim',[-pi pi],'YLim',[0 pi],'ZLim',[-pi pi])
% % set(ax,'XLim',[-pi pi],'YLim',[0 pi],'ZLim',[-pi pi],'DataAspectRatio',[1 1 1])
% box on
% xlabel('phi')
% ylabel('theta')
% zlabel('psi')
%% 
% %% general case
beta2=deg2rad(30);
gamma3=deg2rad(60);
beta3=deg2rad(30);

R=1; L=1;
Lde=simplify(subs(det(Jdkp)));

fx=@(phi,theta,psi)-sin(theta).^2.*(sin(phi).^2 - sin(psi).^2);
interval = [-pi pi 0 pi -pi pi];
fimplicit3(fx,interval,'EdgeColor','none','FaceAlpha',.8)
figure(2)
[phi,theta,psi] = meshgrid(linspace(-pi, pi), ...
                       linspace(0,pi), ...
                       linspace(-pi, pi));
f2 = -sin(theta).^2.*(sin(phi).^2 - sin(psi).^2);
hel = isosurface(phi,theta,psi, f2, 0);

patch(hel,'FaceColor',[.5 1 .5],'EdgeColor','none');
% patch(hel,'FaceColor',[1 .5 0],'EdgeColor','none');
view(3)
camlight
ax = gca;
set(ax,'XLim',[-pi pi],'YLim',[0 pi],'ZLim',[-pi pi],'DataAspectRatio',[1 1 1])
box on
xlabel('phi')
ylabel('theta')
zlabel('psi')
% 
figure(2)
[phi,theta,psi] = meshgrid(linspace(-pi, pi), ...
                        linspace(0.5,pi-0.5), ...
                        linspace(-pi, pi));
f2=((3*3^(1/2).*cos(phi).^2)/8 - (3.*sin(2.*psi))/16 - (3*3^(1/2).*cos(psi).^2)/8 - (3*3^(1/2).*cos(phi).^2.*cos(theta).^2)/8 + (3*3^(1/2).*cos(psi).^2.*cos(theta).^2)/8 + (3.*cos(psi).*cos(theta).^2.*sin(psi))/8 + (3.*cos(psi).*cos(theta).*sin(theta))/8 - (3.*cos(phi).^2.*cos(psi).*cos(theta).*sin(theta))/8 + (3^(1/2).*cos(theta).*sin(psi).*sin(theta))/8 + (3.*cos(phi).*sin(phi).*sin(psi).*sin(theta))/8 + (3^(1/2).*cos(phi).*cos(psi).*sin(phi).*sin(theta))./8 + (3^(1/2).*cos(phi).^2.*cos(theta).*sin(psi).*sin(theta))/8);
%  interval = [-pi pi 0 pi -pi pi];
% fimplicit3(f2,interval,'EdgeColor','none','FaceAlpha',.8)
hel1 = isosurface(phi,theta,psi,f2, 0);
patch(hel1,'FaceColor',[1 .5 0],'EdgeColor','none');

view(3)
camlight
ax = gca;
set(ax,'XLim',[-pi pi],'YLim',[0 pi],'ZLim',[-pi pi])
% set(ax,'XLim',[-pi pi],'YLim',[0 pi],'ZLim',[-pi pi],'DataAspectRatio',[1 1 1])
box on
xlabel('phi')
ylabel('theta')
zlabel('psi')
%% 
% 
% % figure(3)
% % l1=(sqrt((b1-a1)'*(b1-a1)));
% % l2=(sqrt((b2-a2)'*(b2-a2)));
% % l3=(sqrt((b3-a3)'*(b3-a3)));