[x,y]=ode23('f9_10',[0.1:0.05:pi/2],[0.009967,0.99017]);
%poly2str(ployfit(x,y(:,2),3),'t')求微分方程的多项式
x;
yx1=y(:,1);
yl1=y(:,2);
u1=yx1.*yl1.*yl1;
plot(x,u1,'.-b')%图像u(t)
axis([0.1,pi/2,0,0.01]);
xlabel('t'),ylabel('u(t)');
set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]); 
%legend('\alpha=1','真解')
legend('Approximate solution(\alpha=1)')
figure
hold on

plot(x,yx1,'.-r')%图像x(t)
axis([0.1,pi/2,0.001,0.025]);
xlabel('t'),ylabel('x(t)');
set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]); 
%legend('\alpha=1','真解')
legend('Approximate solution(\alpha=1)')
figure
hold on

plot(x,yl1,'.-g')%图像\lambda(t)
%plot(x,y(:,2),'k^',t,y11,'g^',t,y111,'r*')
axis([0.1 pi/2 -1 1]);
xlabel('t'),ylabel('\lambda(t)')
set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]);
%legend('\alpha=1','真解')
legend('Approximate solution(\alpha=1) ')
figure 
hold on

plot(x,u1,'.-b')%图像u(t)
axis([0.1,pi/2,0,0.01]);
xlabel('t'),ylabel('u(t)');
set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]); 
%legend('\alpha=1','真解')
legend('Approximate solution(\alpha=1)')

t=[0.1:0.05:pi/2];
yy=yx1+u1;
trapz(t,yy)
x1=(19735992192100463*cos(t' - atan(990017/1000000) - 1/10).^2)/1000000000000000000;
error=abs(x1-yx1)./yx1;
t=[0.1:0.05:pi/2];
y111=-tan(t' - atan(990017/1000000) - 1/10);
error1=abs(y111-yx1)./yx1;
t=[0.1:0.05:pi/2];
y2=x1.*y111.*y111;
error2=abs(u1-y2)./y2;





