function [dxg,dyg] = dynamicTrackForce(vx,vy,omegaz,omegat,xg,yg)
diffx = abs(x(1,2)-x(1));
diffy = abs(y(2)-y(1));

jx = xg-x;
jy = yg-y;
j = hypot(jy,jx);

vjx = vx - omegaz*y - omegat*r;
vjy = vy + omegaz*x;

delta = atan2(-vjy,-vjx);

dF = sigma*mu*(1-exp(-j/K));
dFx = dF.*cos(delta);
dFy = dF.*sin(delta);
dML = -dFx.*y;
dMr = dFy.*x;

duotrapz = @(d) trapz(xl,trapz(yl,d));

Fx = duotrapz(dFx);
Fy = duotrapz(dFy);
ML = duotrapz(dML);
Mr = duotrapz(dMr);

dxgl = -vx+yg.*omegaz+gradient(xg,diffx)*omegat;
dygl = -vy-xg.*omegaz+gradient(yg,diffy)*omegat;

dxg = [dxgl;dxgr];
dyg = [dygl;dygr];

end