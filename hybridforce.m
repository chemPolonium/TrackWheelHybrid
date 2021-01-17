function hybridforce(vx,vy,omegaz,fx1,fx2)
% HYBRIDFORCE  Get the total force of the track-tire hybrid vehicle.
%   HYBRID(vx,vy,omegaz,fx1,fx2)
%   the hybrid force can calculate the force of the tracks and tires.

[Fxo,Fyo,MLo,Mro,Fxi,Fyi,MLi,Mri] = trackedforce;



    function [Fxr,Fyr,MLr,Mrr,Fxl,Fyl,MLl,Mrl] = trackedforce
        
        % using the coordinate system in Wong's book.
        % where the heading direction is y,
        % the lateral is x,
        % z = cross(x,y) to the sky.
        
        % 履带受力
        % Rx = 10; % R''
        B = 3; % 两条履带中心线之间的距离
        cx = 0; % X方向（车辆横向lateral）上质心到几何中心距离，质心在几何中心左时为正
        cy = 0.1; % Y方向（车辆纵向longitudinal）上质心到几何中心距离，质心在几何中心右时为正
        l = 2; % 履带长度
        % so = 0.5; % Y方向（车辆纵向）上质心到旋转中心距离
        omegaz = 0.1; % 车身旋转速度
        r = 0.3; % 履带轮半径
        omegao = 2; % 外侧履带轮转动速度
        omegai = 1.4; % 内测履带轮转动速度
        b = 0.3; % 履带宽度
        sigmao = 1; % 外侧履带法向压力（压强）
        sigmai = 1; % 内侧履带法向压力（压强）
        mu = 0.5; % 履带和地面的摩擦系数
        K = 0.2;
        
        Rx = vx/omegaz + cx;
        so = vy/omegaz;
        
        x1 = linspace(-b/2,b/2,10);
        y1 = linspace(-l/2-(so-cy),l/2-(so-cy),50)';
        
        x2 = linspace(-b/2,b/2,10);
        y2 = linspace(-l/2-(so-cy),l/2-(so-cy),50)';
        
        jxo = (Rx+B/2+cx+x1).*(cos((l/2+cy-so-y1)*omegaz/r/omegao)-1)-y1.*sin((l/2+cy-so-y1)*omegaz/r/omegao);
        jyo = (Rx+B/2+cx+x1).*sin((l/2+cy-so-y1)*omegaz/r/omegao)-(l/2+cy-so)+y1.*cos((l/2+cy-so-y1)*omegaz/r/omegao);
        jo = hypot(jxo,jyo);
        
        vjyo = (Rx+B/2+cx+x1)*omegaz-r*omegao;
        vjxo = -y1*omegaz;
        vjo = hypot(jxo,jyo);
        
        sindelta1 = vjyo./vjo;
        cosdelta1 = vjxo./vjo;
        
        dFxo = -sigmao*mu*(1-exp(-jo/K)).*cosdelta1;
        dFyo = -sigmao*mu*(1-exp(-jo/K)).*sindelta1;
        dMLo = -(B/2+x1)*sigmao*mu.*(1-exp(-jo/K)).*sindelta1;
        dMro = -y1*sigmao*mu.*(1-exp(-jo/K)).*cosdelta1;
        
        jxi = (Rx-B/2+cx+x2).*(cos((l/2+cy-so-y2)*omegaz/r/omegai)-1)-y2.*sin((l/2+cy-so-y2)*omegaz/r/omegai);
        jyi = (Rx-B/2+cx+x2).*sin((l/2+cy-so-y2)*omegaz/r/omegai)-(l/2+cy-so)+y2.*cos((l/2+cy-so-y2)*omegaz/r/omegai);
        ji = hypot(jxi,jyi);
        
        vjyi = (Rx-B/2+cx+x2)*omegaz-r*omegai;
        vjxi = -y2*omegaz;
        vji = hypot(jxi,jyi);
        
        sindelta2 = vjyi./vji;
        cosdelta2 = vjxi./vji;
        
        dFxi = -sigmai*mu*(1-exp(-ji/K)).*cosdelta2;
        dFyi = -sigmai*mu*(1-exp(-ji/K)).*sindelta2;
        dMLi = -(B/2-x2)*sigmai*mu.*(1-exp(-ji/K)).*sindelta2;
        dMri = -y2*sigmai*mu.*(1-exp(-ji/K)).*cosdelta2;
        
        duotrapz = @(d) trapz(x1,trapz(y1,d));
        
        Fxr = duotrapz(dFxo);
        Fyr = duotrapz(dFyo);
        MLr = duotrapz(dMLo);
        Mrr = duotrapz(dMro);
        
        Fxl = duotrapz(dFxi);
        Fyl = duotrapz(dFyi);
        MLl = duotrapz(dMLi);
        Mrl = duotrapz(dMri);
        
        % 需要使用内置积分器时，创建这样一个函数
        % dFxo = getdFxo(x1,y1);
        
        % 采用matlab自带的积分器，速度相对慢一点
        % 需要把中间一部分运算创建成一个函数
        % Fxo = integral2(@(x,y) getdFxo(x,y),-b/2,b/2,-l/2-(so-cy),l/2-(so-cy))
        
    end

    function Fy = fialatireforce(talpha,Fx)
        % 轮胎受力
        % 采用附着椭圆，fiala模型
        % tire forces
        % using fiala model, which uses a 3rd order tylor expansion to fit the
        % force of a tire.
        mu = 0.55;
        Ca = 9e4;
        Fz = 1e3;
        eta = sqrt(mu^2*Fz^2-Fx^2)/mu/Fz;
        tasl = 3*eta*mu*Fz/Ca;
        if abs(talpha) < tasl
            Fy = -Ca*talpha+...
                Ca^2/(3*eta*mu*Fz)*abs(talpha)*talpha-...
                Ca^3/(27*eta^2*mu^2*Fz^2)*talpha^3;
        else
            Fy = -eta*mu*Fz*sign(talpha);
        end
    end

    function [Fx,Fy,Mz] = tireforce(l,delta,FxTire)
        % delta: steer angle
        talpha = (vy+omegaz*l)/vx;
        FyTire = tireforce(talpha,FxTire);
        Fx = FxTire*cos(delta)-FyTire*sin(delta);
        Fy = FyTire*cos(delta)+FxTire*sin(delta);
        Mz = cx*Fx+l*Fy;
    end

end
