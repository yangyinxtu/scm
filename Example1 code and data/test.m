clear; clc;
format long e
A=0;
B=1;
%alpha=0.999999999;
%alpha=0.9;
%alpha=0.8;
%alpha=0.7;
alpha=0.999999;
for N = 2 : 2 :20
   
    [theta1, w1] = jacobi_compute(N+1,0.0,0.0,-1,1);
    theta1 = theta1';  w1 = w1';
    
    [theta2, w2] = jacobi_compute(N+1,-alpha,0.0,-1,1);
    
    theta2 = theta2';
    x=theta2; 
    w2 = w2';
    
    [theta3, w3] = jacobi_compute(N+1,0.0,0.0,-1,1);
    theta3 = theta3'; w3 = w3';
    
    [theta4, w4] = jacobi_compute(N+1,0.0,-alpha,-1,1);
    theta4 = theta4';  theta=theta4; w4 = w4';
   
    t=0.5*x+1/2;
    a = (1+x)/2; b = (x-1)/2;
    c = (1-x)/2; d = (x+1)/2;
    
    for i=1:N+1
        for j=1:N+1
            bb=0;
            %计算积分
            for k=1:N+1
                %F_1(x,j,N,a(i)*theta(k)+b(i))计算第j个配置点上的N次插值多项式
                %在a(i)*theta(k)+b(i)处的值�?实现见F_1.m
                
                bb=bb + w1(k) *F_1(x,j,N,a(i)*theta1(k)+b(i));
            end
            %
            M(i,j)=1/2*((x(i)+1)/2) * bb;
        end
    end


    %计算矩阵L
    for i=1:N+1
        for j=1:N+1
            ll=0;
            %计算积分
            for k=1:N+1
                %F_1(x,j,N,a(i)*theta(k)+b(i))计算第j个配置点上的N次插值多项式
                %在a(i)*theta(k)+b(i)处的值�?实现见F_1.m
                %kxs=1;
                 %kxs=u((a(i)*theta2(k)+b(i)));
                ll=ll + w2(k) *F_1(x,j,N,a(i)*theta2(k)+b(i));
            end
            %L(i,j)= -(B/(gamm(0.5))*(1/2)^(0.5)*((x(i)+1)/2)^(1-alpha) * ll);
            L(i,j)= B/(gamma(1-alpha))*(1/2).^(1-alpha)*((x(i)+1)/2).^(1-alpha) * ll;
        end
    end
     for i=1:N+1
        for j=1:N+1
            cc=0;
            %计算积分
            for k=1:N+1
                %F_1(x,j,N,a(i)*theta(k)+b(i))计算第j个配置点上的N次插值多项式
                %在a(i)*theta(k)+b(i)处的值�?实现见F_1.m
                
                cc=cc + w3(k) *F_1(x,j,N,c(i)*theta3(k)+d(i));
            end
            %
            Q(i,j)= 1/2*((1-x(i))/2) * cc;
        end
    end


    %计算矩阵L
    for i=1:N+1
        for j=1:N+1
            dd=0;
            %计算积分
            for k=1:N+1
                %F_1(x,j,N,a(i)*theta(k)+b(i))计算第j个配置点上的N次插值多项式
                %在a(i)*theta(k)+b(i)处的值�?实现见F_1.m
                %kxs=1;
                dd=dd + w4(k)*F_1(x,j,N,c(i)*theta4(k)+d(i));
            end
            %R(i,j)= -(B/(gamm(1-alpha))*(1/2)^(1-alpha).*((1-x(i))/2)^(1-alpha) * dd);
            R(i,j)= B/(gamma(1-alpha))*(1/2)^(1-alpha).*((1-x(i))/2)^(1-alpha) * dd;
        end
    end
     
 %K=A+R+Q;
 %J=A+L+M;
% KJ=[K,M;
%         J,-Q];
%u_1=ones(N+1,1);
%u_1=-1*u_1; FF=[u_1;u_1]; uz=KJ\FF;dlamda=[uz(1:N+1)]; du=[uz(N+2:2*(N+1))];lamda=-Q*dlamda; u=M*du+1;

K=A+R;
J=A+L;
%for i=1:N+1
 %   ab(i)=1;
%end
%E=diag(ab);
E=eye(N+1);

O=zeros(N+1);
KJ=[K,-E,E,O;
    O,E,E,J;
    O,E,O,-M;
    -Q,O,E,O];
u_1=ones(N+1,1);
D=zeros(N+1,1);
FF=[D;D;u_1;D];
uz=KJ\FF;
dlamda=[uz(1:N+1)];
u=[uz(N+2:2*(N+1))];
lamda=[uz(2*N+3:3*(N+1))];
du=[uz(3*N+4:4*(N+1))];
f=-lamda;
df=-dlamda;
beta=-(cosh(sqrt(2))+sqrt(2)*sinh(sqrt(2)))/(sqrt(2)*cosh(sqrt(2))+sinh(sqrt(2)));
%  t
 index=find(t>=0);
 t=t(index);
 exact= cosh(sqrt(2).*t)+beta*sinh(sqrt(2).*t);
 
 H1=0:0.01:1;
 u_h=zeros(size(H1'));
 for j=1:N+1
     u_h=u_h+u(j)*F_1(t,j,N,H1');
 end
 exact_h=cosh(sqrt(2).*H1')+beta*sinh(sqrt(2).*H1');
%  AAA=size(exact_h)
%  BBB=size(H')
 
error(N/2) = max( abs( u_h - exact_h) );
L_Norm = w2'.*(u - exact').^2;
L_error(N/2)=sqrt(sum(L_Norm));

index=find(t>=0);
 t=t(index);
fexact=((1+sqrt(2)*beta)*cosh(sqrt(2)*t)+(sqrt(2)+beta)*sinh(sqrt(2)*t));
H1=0:0.01:1;
f_h=zeros(size(H1'));
for j=1:N+1
    f_h=f_h+f(j)*F_1(t,j,N,H1');
end
fexact_h=((1+sqrt(2)*beta)*cosh(sqrt(2)*H1')+(sqrt(2)+beta)*sinh(sqrt(2)*H1'));

ferror(N/2)=max(abs(f_h-fexact_h));
fL_Norm=w2'.*(f-fexact').^2;
fL_error(N/2)=sqrt(sum(fL_Norm));


index=find(t>=0);
 t=t(index);
lamdaexact=-((1+sqrt(2)*beta)*cosh(sqrt(2)*t)+(sqrt(2)+beta)*sinh(sqrt(2)*t));
H1=0:0.01:1;
lamda_h=zeros(size(H1'));
for j=1:N+1
    lamda_h=lamda_h+lamda(j)*F_1(t,j,N,H1');
end
lamdaexact_h=-((1+sqrt(2)*beta)*cosh(sqrt(2)*H1')+(sqrt(2)+beta)*sinh(sqrt(2)*H1'));

lamdaerror(N/2)=max(abs(lamda_h-lamdaexact_h));
lamdaL_Norm=w2'.*(lamda-lamdaexact').^2;
lamdaL_error(N/2)=sqrt(sum(lamdaL_Norm));
%����\lamda(t)�ĵ������
index=find(t>=0);
 t=t(index);
 %- (6796620497027857*2^(1/2)*cosh(2^(1/2)*t))/70368744177664 - (4841120814688713*2^(1/2)*sinh(2^(1/2)*t))/35184372088832
dlamdaexact=-(1955876548596813*2^(1/2)*cosh(2^(1/2)*t))/4503599627370496 - (434393121504351*2^(1/2)*sinh(2^(1/2)*t))/1125899906842624;
H1=0:0.01:1;
dlamda_h=zeros(size(H1'));
for j=1:N+1
    dlamda_h=dlamda_h+lamda(j)*F_1(t,j,N,H1');
end
dlamdaexact_h=-(1955876548596813*2^(1/2)*cosh(2^(1/2)*H1'))/4503599627370496 - (434393121504351*2^(1/2)*sinh(2^(1/2)*H1'))/1125899906842624;

dlamdaerror(N/2)=max(abs(dlamda_h-dlamdaexact_h));
dlamdaL_Norm=w2'.*(lamda-dlamdaexact').^2;
dlamdaL_error(N/2)=sqrt(sum(dlamdaL_Norm));


index=find(t>=0);
 t=t(index);
 %- (6796620497027857*2^(1/2)*cosh(2^(1/2)*t))/70368744177664 - (4841120814688713*2^(1/2)*sinh(2^(1/2)*t))/35184372088832
dfexact=(1955876548596813*2^(1/2)*cosh(2^(1/2)*t))/4503599627370496 - (434393121504351*2^(1/2)*sinh(2^(1/2)*t))/1125899906842624;
H1=0:0.01:1;
df_h=zeros(size(H1'));
for j=1:N+1
    df_h=df_h+f(j)*F_1(t,j,N,H1');
end
dfexact_h=(1955876548596813*2^(1/2)*cosh(2^(1/2)*H1'))/4503599627370496 - (434393121504351*2^(1/2)*sinh(2^(1/2)*H1'))/1125899906842624;

dferror(N/2)=max(abs(df_h-dfexact_h));
dfL_Norm=w2'.*(f-dfexact').^2;
dfL_error(N/2)=sqrt(sum(dfL_Norm));

index=find(t>=0);
 t=t(index);

 dexact=2^(1/2)*sinh(2^(1/2)*t) - (34477930655695*2^(1/2)*cosh(2^(1/2)*t))/35184372088832;
 H1=0:0.01:1;
 du_h=zeros(size(H1'));
 for j=1:N+1
     du_h=du_h+du(j)*F_1(t,j,N,H1');
 end
 dexact_h=2^(1/2)*sinh(2^(1/2)*H1') - (34477930655695*2^(1/2)*cosh(2^(1/2)*H1'))/35184372088832;
 
    derror(N/2) = max( abs( du_h - dexact_h ) );
    dL_Norm = w4'.*(du - dexact').^2;
    dL_error(N/2) =  sqrt(sum(dL_Norm) );    
end  
yy=1/2*(u.^2+f.^2);
trapz(t,yy)
%下面是将结果用图来显示，分别画出真解，数值解及误差图�?
plot(t,u,'b*',H1',exact_h,'m','MarkerSize',10)
%plot(t,u,'k:','MarkerSize',10)
%plot(t,u,'g^-','MarkerSize',10)
%plot(t,u,'k--','MarkerSize',10)
hold on
xlabel('t'),ylabel('x(t)')
figure
plot(t,lamda,'b*',H1,lamdaexact_h,'m','MarkerSize',10)
hold on
xlabel('t'),ylabel('\lambda(t)')
legend('\alpha=1','���')
figure
plot(t,f,'b*',H1,fexact_h,'m','MarkerSize',10)
 hold on
 xlabel('t'),ylabel('u(t)')
 %legend('Approximate solution  ','Exact solution  ')
 legend('\alpha=1','\alpha=0.99')
%legend('Approximate solution  ','Exact solution  ')
 figure
%plot(t,du,'g*',H1,dexact_h,'r','MarkerSize',10)
%hold on

%legend('du Approximate Solution','du Exact Solution',2)
%legend('Approximate derivative  ','Exact derivative  ')
%figure 

semilogy(2:2:20,error,'bo-.','LineWidth',2,'MarkerSize',10)
hold on
%semilogy(2:2:20,L_error,'mdiamond-.','LineWidth',2,'MarkerSize',10)
%semilogy(2:2:20,dL_error,'rp--','LineWidth',2,'MarkerSize',10)
semilogy(2:2:20,derror,'g^-.','LineWidth',2,'MarkerSize',10)
axis tight;
xlabel('t(2\leq N\leq 20)'),ylabel('x(t)-x^{N}(t)')
legend('L^{\infty}-error  ','L^{2}_{\omega^{-\mu,-\mu}}-error ')

figure
%xlabel('$\dot{x_j}(t)$','Interpreter','latex')


%%semilogy(2:2:20,error,'bo-.','LineWidth',2,'MarkerSize',10)
%semilogy(2:2:20,L_error,'mdiamond-.','LineWidth',2,'MarkerSize',10)
%%semilogy(2:2:20,derror,'g^-.','LineWidth',2,'MarkerSize',10)
%hold on
%semilogy(2:2:20,dL_error,'rp--','LineWidth',2,'MarkerSize',10)
%axis tight;
%xlabel('t(2\leq N\leq 20)'),ylabel('$\dot{x}(t)$-$\dot{x}^{N}(t)$','Interpreter','latex')
%legend('L^{\infty}-error  ','L^{2}_{\omega^{-\mu,-\mu}}-error  ')

%figure
semilogy(2:2:20,ferror,'g^-.','LineWidth',2,'MarkerSize',10)
hold on 
%semilogy(2:2:20,dfL_error,'rp--','LineWidth',2,'MarkerSize',10)
%semilogy(2:2:20,dferror,'bo-.','LineWidth',2,'MarkerSize',10)
semilogy(2:2:20,fL_error,'mdiamond-.','LineWidth',2,'MarkerSize',10)
axis tight;

xlabel('t(2\leq N\leq 20)'),ylabel('$\lambda (t)-\lambda^{N}(t)$','Interpreter','latex')
legend('L^{\infty}-error  ','L^{2}_{\omega^{-\mu,-\mu}}-error  ')


%ylabel('||u-u^N||_{\infty}')
%title('','Fontsize',10)
set(gca,'FontSize',10)
%print -depsc 1-D.eps % save figure to eps file
                
 
