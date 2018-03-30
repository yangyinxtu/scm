clc;clear;
a=0.1;b=pi/2;m=30;
x=cos(pi*(0:m)'/m);
x=(b+a)/2+(b-a)*x/2;
Dx=barylag(x,3);
Dx1=Dx(:,:,1);
Dx2=Dx(:,:,2);
Dx3=Dx(:,:,3);
m=length(x);
I=eye(m);
II=eye(2*m);
%u=inline('1/3/pi*sin(3*pi*x)','x');
%v=inline('3*pi*x-3*pi/2','x');
e=10e-10;
u0=zeros(m,1);
v0=zeros(m,1);
IE=zeros(m,1);
U0=[u0;v0];
%L=[Dx1,0*I;0*I,Dx1];
%F=[-(1+u0.^2);2*u0.*v0];
L=[Dx1+2*diag(u0),0*I;0*I,Dx1-2*diag(u0)];
F=[-1+[u0.^2];IE];
L(1,:)=II(1,:);
F(1)=0.009769;
L(m+1,:)=II(m+1,:);
F(m+1)=0.009967;
U1=L\F;
mm=1;
while norm(U1-U0,inf)>e
    U0=U1;
    mm=mm+1;
    u0=U0(1:m);
    v0=U0(m+1:end);
    %L=[Dx1,0*I;0*I,Dx1];
    %F=[-(1+u0.^2);2*u0.*v0];
    L=[Dx1+2*diag(u0),0*I;0*I,Dx1-2*diag(u0)];
    F=[-1+[u0.^2];IE];
    L(1,:)=II(1,:);
    F(1)=0.009769;
    L(m+1,:)=II(m+1,:);
    F(m+1)=0.009967;
    U1=L\F;
    if mm>100 break;
    end
end
u1=U1(1:m)
v1=U1(m+1:end)
n1=v1.*u1.^2
%ue=u(x);
%ve=v(x);
%err=[m mm,norm(u1-ue),norm(u1-ue)/norm(ue)...
   % norm(v1-ve),norm(v1-ve)/norm(ve)]
   

plot(x,u1,'.-k');
%axis([0.1 pi/2 0 10 ]);
%xlabel('t'),ylabel('\lambda(t)');
%set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]); 
%legend('numerical solution')
figure
hold on;

plot(x,v1,'.-k');
%axis([0.1 pi/2 0 0.01]);
%xlabel('t'),ylabel('x(t)');
%set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]); 
%legend('numerical solution ')
figure
hold on



plot(x,n1,'.-b')
axis([0.1 pi/2 0 0.01]);
xlabel('t'),ylabel('u(t)');
set(gca, 'XTick', [0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 pi/2]); 
legend('Approximate solution  ')




err=[Dx1*u1+(1+u0.^2)...
    Dx1*v1-2*u0.*v0]
%plot(x,err,'+')