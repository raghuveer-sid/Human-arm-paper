% ----- Function Subroutine Below ----- %
function [x,y,z] = torusmod(r,d,p)
[th,phi] = meshgrid(-3.1:0.1:3.2);
x = (r*cos(th/p)+d).*cos(phi);
y = (r*cos(th/p)+d).*sin(phi);
z = r/2*sin(th/p);
end

%%