function ff= F_1(x,j,N,s)
%LAGRANGE_INTERPOLATION_POLY Summary of this function goes here
%   Detailed explanation goes here
ff=ones(size(s));
for m=1:N+1
    if j==m
        ff=ff;
    else ff=ff.*(s-x(m))/(x(j)-x(m));
    end

end

