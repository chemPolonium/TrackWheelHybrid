% virtual y accelaration
% virtual z angle velocity
ayvlist = 2:0.1:2.5;
avzvlist = 0.1:0.1:0.6;
[ayv,avzv] = meshgrid(ayvlist,avzvlist);

vsl = (sqrt(ayv./avzv)-sqrt(ayv.*avzv))/2;
vsr = (sqrt(ayv./avzv)+sqrt(ayv.*avzv))/2;

% x velocity, y velocity, z angle velocity
[vx,vy,avz] = deal(zeros(size(vsl)));

set_param("trkwhlmdldynf/XY Graph","Commented","on");
for i = 1:numel(vsl)
    if vsl(i) == vsr(i)
        vx(i) = vsl(i);
        vy(i) = 0;
        avz(i) = 0;
        continue
    end
    ivsl = vsl(i);
    ivsr = vsr(i);
    disp("sim for vsl = "+ivsl+", vsr = "+ivsr);
    out = sim("trkwhlmdldynf.slx");
    ivx = out.vel.vx.Data(end);
    ivy = out.vel.vy.Data(end);
    iomegaz = out.vel.omegaz.Data(end);
    vx(i) = ivx;
    vy(i) = ivy;
    avz(i) = iomegaz;
end
set_param("trkwhlmdldynf/XY Graph","Commented","off");

yc = vx./avz;
yl = (vx-vsl)./avz;
yr = (vx-vsr)./avz;
xc = -vy./avz;

figurePosition = [100,500,500,400];
for varName = ["xc" "yl" "yr"]
figure("Position",figurePosition);
surf(ayv,avzv,eval(varName));
xlabel("ayv");
ylabel("avzv");
title(varName);
exportgraphics(gca,"pic\dynf_"+varName+".png");
figurePosition(1) = figurePosition(1)+figurePosition(3);
end
