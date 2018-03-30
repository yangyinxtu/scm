clear; clc;
format long e
A=0;
B=1;
alpha=0.9999999999;
%alpha=0.9;
%alpha=0.8;
%alpha=0.7;

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
            L(i,j)= B/(gamma(1-alpha))*(1/2)^(1-alpha)*((x(i)+1)/2)^(1-alpha) * ll;
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
     
 
K=A+R;
J=A+L;
E=eye(N+1);
O=zeros(N+1);
KJ=[K,-E,E,O;
    O,E,E,J;
    O,E,O,-M;
    -Q,O,E,O];
u_1=10*ones(N+1,1);
D=zeros(N+1,1);
FF=[D;D;u_1;D];
uz=KJ\FF;
dlamda=[uz(1:N+1)];
u=[uz(N+2:2*(N+1))];
lamda=[uz(2*N+3:3*(N+1))];
du=[uz(3*N+4:4*(N+1))];
f=-lamda;

  index=find(t>=0);
  t=t(index);
  exact= 0.1006126842*exp(sqrt(2).*t)+9.899387312*exp(-sqrt(2).*t);
 
  H1=0:0.01:1;
  u_h=zeros(size(H1'));
  for j=1:N+1
     u_h=u_h+u(j)*F_1(t,j,N,H1');
  end
  exact_h=0.1006126842*exp(sqrt(2).*H1')+9.899387312*exp(-sqrt(2).*H1');

  error(N/2) = max( abs( u_h - exact_h) );
  L_Norm = w2'.*(u - exact').^2;
  L_error(N/2)=sqrt(sum(L_Norm));


  index=find(t>=0);
  t=t(index);

  dexact=(3624953976574923*2^(1/2)*exp(2^(1/2)*t))/36028797018963968 - (5572859626189927*2^(1/2)*exp(-2^(1/2)*t))/562949953421312;

  H1=0:0.01:1;
  du_h=zeros(size(H1'));
  for j=1:N+1
     du_h=du_h+du(j)*F_1(t,j,N,H1');
  end
  dexact_h=(3624953976574923*2^(1/2)*exp(2^(1/2)*H1'))/36028797018963968 - (5572859626189927*2^(1/2)*exp(-2^(1/2)*H1'))/562949953421312;
 
    derror(N/2) = max( abs( du_h - dexact_h ) );
    dL_Norm = w4'.*(du - dexact').^2;
    dL_error(N/2) =  sqrt(sum(dL_Norm) );
 
  index=find(t>=0);
  t=t(index);
  fexact=0.1006126842*(sqrt(2)+1)*exp(sqrt(2).*t)-9.899387312*(sqrt(2)-1)*exp(-sqrt(2).*t);
  H1=0:0.01:1;
  f_h=zeros(size(H1'));
  for j=1:N+1
    f_h=f_h+f(j)*F_1(t,j,N,H1');
  end
  fexact_h=0.1006126842*(sqrt(2)+1)*exp(sqrt(2).*H1')-9.899387312*(sqrt(2)-1)*exp(-sqrt(2).*H1');

  ferror(N/2)=max(abs(f_h-fexact_h));
  fL_Norm=w2'.*(f-fexact').^2;
  fL_error(N/2)=sqrt(sum(fL_Norm));

  index=find(t>=0);
  t=t(index);
  lamdaexact=(-0.1006126842*(sqrt(2)+1)*exp(sqrt(2).*t)+9.899387312*(sqrt(2)-1)*exp(-sqrt(2).*t));
  H1=0:0.01:1;
  lamda_h=zeros(size(H1'));
  for j=1:N+1
    lamda_h=lamda_h+lamda(j)*F_1(t,j,N,H1');
  end
  lamdaexact_h=(-0.1006126842*(sqrt(2)+1)*exp(sqrt(2).*H1')+9.899387312*(sqrt(2)-1)*exp(-sqrt(2).*H1'));

  lamdaerror(N/2)=max(abs(lamda_h-lamdaexact_h));
  lamdaL_Norm=w1'.*(lamda-lamdaexact').^2;
  lamdaL_error(N/2)=sqrt(sum(lamdaL_Norm));
  %����\lamda(t)�ĵ������
  index=find(t>=0);
  t=t(index);
  %- (6796620497027857*2^(1/2)*cosh(2^(1/2)*t))/70368744177664 - (4841120814688713*2^(1/2)*sinh(2^(1/2)*t))/35184372088832
  dlamdaexact=(-8751413053225461*2^(1/2)*exp(2^(1/2)*t))/36028797018963968 +(2308354038369325*2^(1/2)*exp(-2^(1/2)*t))/562949953421312;
  H1=0:0.01:1;
  dlamda_h=zeros(size(H1'));
  for j=1:N+1
    dlamda_h=dlamda_h+lamda(j)*F_1(t,j,N,H1');
  end
  dlamdaexact_h=(-8751413053225461*2^(1/2)*exp(2^(1/2)*H1'))/36028797018963968 +(2308354038369325*2^(1/2)*exp(-2^(1/2)*H1'))/562949953421312;

  dlamdaerror(N/2)=max(abs(dlamda_h-dlamdaexact_h));
  dlamdaL_Norm=w3'.*(lamda-dlamdaexact').^2;
  dlamdaL_error(N/2)=sqrt(sum(dlamdaL_Norm));   
    
    
    
end  
 plot(t,u,'b*',H1',exact_h,'m','MarkerSize',10)
 hold on
 xlabel('t'),ylabel('x(t)')  
 figure
 
 plot(t,du,'b*',H1',dexact_h,'m','MarkerSize',10)
 hold on
 xlabel('t'),ylabel('$\dot{x}(t)$','Interpreter','latex')
 figure
 
 plot(t,f,'b*',H1,fexact_h,'m','MarkerSize',10)
 hold on
 xlabel('t'),ylabel('u(t)')
 figure
 
 plot(t,lamda,'b*',H1,lamdaexact_h,'m','MarkerSize',10)
 hold on 
 xlabel('t'),ylabel('\lambda(t)')
figure 

semilogy(2:2:20,error,'bo-.','LineWidth',2,'MarkerSize',10)
hold on
semilogy(2:2:20,derror,'g^-.','LineWidth',2,'MarkerSize',10)
axis tight;
xlabel('t(2\leq N\leq 20)'),ylabel('x(t)-x^{N}(t)')
legend('L^{\infty}-error  ','L^{2}_{\omega^{-\mu,-\mu}}-error ')
figure 

semilogy(2:2:20,lamdaerror,'bo-.','LineWidth',2,'MarkerSize',10)
hold on
semilogy(2:2:20,lamdaL_error,'g^-.','LineWidth',2,'MarkerSize',10)
axis tight;
xlabel('t(2\leq N\leq 20)'),ylabel('$\lambda (t)-\lambda^{N}(t)$','Interpreter','latex')
legend('L^{\infty}-error  ','L^{2}_{\omega^{-\mu,-\mu}}-error ')
%ylabel('||u-u^N||_{\infty}')
%title('','Fontsize',10)
set(gca,'FontSize',10)
%print -depsc 1-D.eps % save figure to eps file
                
%g=inline('0.5*u.^2+0.5*f.^2');
%JI=quadl(g,0,1)
 
