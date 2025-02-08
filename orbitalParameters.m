clear; clc;

% find h,e,i,oBig,oSmall,theta
mu = 398600.4418; %km^3s^-2
r = [-7242.74, -1298.97, -1816.09];
v = [1.950, -6.834, 1.405];

dist = sqrt((r(1))^2 + (r(2))^2 + (r(3))^2);
spd = sqrt((v(1))^2 + (v(2))^2 + (v(3))^2);
vr = dot(r,(v/dist));

h = cross(r,v);
hMag = sqrt(dot(h,h)) %h

i = acosd(h(3)/hMag) %i

n = [-h(2), h(1), 0];
if (n(2) < 0)
   oBig = 360 - acosd(n(1)/sqrt(dot(n,n)))
else
    oBig = acosd(n(1)/sqrt(dot(n,n)))
end %Omega

eMag = sqrt(1 + (((hMag^2)/(mu^2)) * ((spd^2)-(2*mu/dist)))) %e

e = (1/mu) *( (((spd^2)-(mu/dist))*r) - (dist*vr*v));

% sanity check -> eMagCheck = sqrt(dot(e,e));

if (e(3) < 0)
    oSmall = 360 - acosd(dot((n/sqrt(dot(n,n))),(e/eMag)))
else
    oSmall = acosd(dot((n/sqrt(dot(n,n))),(e/eMag)))
end %omega

if (vr < 0)
    theta = 360 - acosd(dot((e/eMag),(r/dist)))
else
    theta = acosd(dot((e/eMag),(r/dist)))
end %theta



 
