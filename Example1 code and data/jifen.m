   function z=jifen(x)
h=0.000001;
t=(1+x)/2;
a=h:h:t
u=exp(a)/sqrt(a);
u0=exp(h)/sqrt(h);
un=exp(t)./sqrt(t);
ji=sum(u,2)-u0/2-un/2 ;
z=ji*h;  