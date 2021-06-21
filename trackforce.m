function [Fx,Fy,Mz] = trackforce(vx,vy,omegaz,omegal,omegar)
% HYBRIDFORCE  Get the total force of the track-tire hybrid vehicle.
%   HYBRID(vx,vy,omegaz,fx1,fx2)
%   the hybrid force can calculate the force of the tracks and tires.

% omegaz 车身旋转速度
% omegar 右侧履带轮转动速度
% omegal 左测履带轮转动速度

m = 3e3; % 车身质量

r = 0.3; % 履带轮半径
b = 0.3; % 履带宽度
mu = 0.5; % 履带和地面的摩擦系数
K = 0.2; % 土壤单位长度

B = 3; % 两条履带中心线之间的距离
l = 4; % 履带长度

sigmar = m/(2*l*b); % 右侧履带法向压力（压强）
sigmal = m/(2*l*b); % 左侧履带法向压力（压强）

% xll: x left list
xll = linspace(-l/2,l/2,50);
yll = linspace(-b/2+B/2,b/2+B/2,10)';

xrl = linspace(-l/2,l/2,50);
yrl = linspace(-b/2-B/2,b/2-B/2,10)';

[Fxl,Fyl,MLl,Mrl] = singleTrackForce(xll,yll,omegal,sigmal);
[Fxr,Fyr,MLr,Mrr] = singleTrackForce(xrl,yrl,omegar,sigmar);

Fx = Fxl+Fxr;
Fy = Fyl+Fyr;
Mz = MLl+Mrl+MLr+Mrr;

    function [Fx,Fy,ML,Mr] = singleTrackForce(xl,yl,omega,sigma)
        [x,y] = meshgrid(xl,yl);
        [jx,jy] = vel2dis(-vx,-vy,-omegaz,(l/2-x)./omega,l/2,y);
        jx = jx - x;
        jy = jy - y;
        j = hypot(jy,jx);
        
        vjx = vx - omegaz*y - omega*r;
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
    end

end
