function [H,x]=createH(n)
V0 = 1;
delta = 3;
a = 50;
dx = a/n;
opts.dx = dx;
x = [0:dx:a-dx]';
V_kp = zeros(n,1);
V_free = V_kp;
for j = 0:5
    x0 = x - a*j/(5);
    tempV = -V0*exp(-(x0).^2/(delta^2));
    V_kp = V_kp + tempV;
end

l1=n;

Lap_1D = -2*speye(l1,l1) + spdiags(ones(l1-1,1),-1,l1,l1) + spdiags(ones(l1,1),1,l1,l1);
Lap_1D(1,l1) = 1;  Lap_1D(l1,1) = 1;
H = -1/2*Lap_1D/dx/dx;% + spdiags(V_kp,0,l1,l1);
% H = -1/2*Lap_1D/dx/dx + spdiags(V_kp,0,l1,l1);
end