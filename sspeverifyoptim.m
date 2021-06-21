omegallist = 19:0.2:21;
omegarlist = 19:0.2:21;
[omegal,omegar] = meshgrid(omegallist,omegarlist);

[vx,vy,omegaz] = arrayfun(@solveeq2,omegal,omegar);

function [vx,vy,omegaz] = solveeq(omegal,omegar)
[x,fval] = fsolve(@(x) eqfcn(x(1),x(2),x(3),omegal,omegar),[4 0 0]);
vx = x(1);
vy = x(2);
omegaz = x(3);
disp(norm(fval));
end

function [vx,vy,omegaz] = solveeq2(omegal,omegar)
[x,fval] = fminsearch(@(x) eqfcn2(x(1),x(2),x(3),omegal,omegar),[4 0 0]);
vx = x(1);
vy = x(2);
omegaz = x(3);
if fval > 0.1
    disp("not solved. omegal:"+string(omegal)+", omegar:"+string(omegar)+...
        " x:"+string(x(1))+","+string(x(2))+","+string(x(3))+", fval:"+string(fval));
end
end

function c = eqfcn(vx,vy,omegaz,omegal,omegar)
    [Fx,Fy,Mz] = trackforce(vx,vy,omegaz,omegal,omegar);
    m = 3e3;
    c = [m.*vx.*omegaz - Fy;m.*vy.*omegaz + Fx;Mz];
end

function c = eqfcn2(vx,vy,omegaz,omegal,omegar)
    [Fx,Fy,Mz] = trackforce(vx,vy,omegaz,omegal,omegar);
    m = 3e3;
    c = (m.*vx.*omegaz - Fy).^2 + (m.*vy.*omegaz + Fx).^2 + Mz.^2;
end