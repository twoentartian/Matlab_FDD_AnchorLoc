figure;
hold on;
T=[-1;0]; R=[1;0]; A=[sqrt(3);0.5];
plot(-1,0,'rp','Markersize',10,'linewidth',2);text(-1,0,'发射机T');
plot(1,0,'gp','Markersize',10,'linewidth',2);text(1,0,'接收机R');
theta=0:pi/100:2*pi; plot(2*cos(theta),sin(theta),'b');
plot(sqrt(3),0.5,'y^','Markersize',10,'linewidth',2); text(sqrt(3),0.5,'转发机A');
plot([-1,sqrt(3),1,-1],[0,0.5,0,0],'black');
axis equal;