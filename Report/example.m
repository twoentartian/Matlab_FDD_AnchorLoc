figure;
hold on;
T=[-1;0]; R=[1;0]; A=[sqrt(3);0.5];
plot(-1,0,'rp','Markersize',10,'linewidth',2);text(-1,0,'发射机T');
plot(1,0,'gp','Markersize',10,'linewidth',2);text(1,0,'接收机R');
theta=0:pi/100:2*pi; plot(2*cos(theta),sin(theta),'b');
plot(sqrt(3),0.5,'y^','Markersize',10,'linewidth',2); text(sqrt(3),0.5,'转发机A');
plot([-1,sqrt(3),1],[0,0.5,0],'b');
%
plot(2,-1,'co','Markersize',10,'linewidth',2);text(2,-1,'监听机S1');
ezplot('sqrt((x+1)^2+y^2)+sqrt((x-2)^2+(y+1)^2)-4.3012',[-3,3,-3,3]);
plot([-1,sqrt(3),2],[0,0.5,-1],'g');
%
plot(1,1,'co','Markersize',10,'linewidth',2);text(1,1,'监听机S2');
ezplot('sqrt((x+1)^2+y^2)+sqrt((x-1)^2+(y-1)^2)-3.6639',[-3,3,-3,3]);
plot([-1,sqrt(3),1],[0,0.5,1]);

axis equal;