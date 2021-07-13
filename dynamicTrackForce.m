function [Fx,Fy,Mz,dxg,dyg] = dynamicTrackForce(vx,vy,omegaz,omegat,sigma,r,mu,K,x,y,xg,yg)
diffx = abs(x(1,2)-x(1));

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

duotrapz = @(d) trapz(x(1,:),trapz(y(:,1),d));

Fx = duotrapz(dFx);
Fy = duotrapz(dFy);
ML = duotrapz(dML);
Mr = duotrapz(dMr);
Mz = ML+Mr;

dxg = -vx+yg.*omegaz+gradient(xg,diffx).*r.*omegat;
dyg = -vy-xg.*omegaz+gradient(yg,diffx).*r.*omegat;

dxg(:,end) = 0;
dyg(:,end) = 0;
