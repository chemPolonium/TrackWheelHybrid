vsllist = 1.5:0.3:3;
vsrlist = 1.5:0.3:3;
[vsl,vsr] = meshgrid(vsllist,vsrlist);

[vx,vy,omegaz] = deal(zeros(size(vsl)));

set_param("trkwhlmdldynf/XY Graph","Commented","on");
for i = 1:numel(vsl)
    if vsl(i) == vsr(i)
        vx(i) = vsl(i);
        vy(i) = 0;
        omegaz(i) = 0;
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
    omegaz(i) = iomegaz;
end
set_param("trkwhlmdldynf/XY Graph","Commented","off");

yc = vx./omegaz;
yl = (vx-vsl)./omegaz;
yr = (vx-vsr)./omegaz;
xc = -vy./omegaz;

abs((vsl-vsr).*(vsl+vsr));
abs((vsl-vsr)./(vsl+vsr));

figurePosition = [100,500,500,400];
for varName = ["vx" "vy" "omegaz"]
figure("Position",figurePosition);
surf(vsl,vsr,eval(varName));
xlabel("vsl");
ylabel("vsr");
title(varName);
exportgraphics(gca,"pic\dynf_"+varName+".png");
figurePosition(1) = figurePosition(1)+figurePosition(3);
end

