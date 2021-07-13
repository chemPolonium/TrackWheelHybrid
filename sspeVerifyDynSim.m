omegallist = 5:1:10;
omegarlist = 5:1:10;
[omegal,omegar] = meshgrid(omegallist,omegarlist);

[vx,vy,omegaz] = deal(zeros(size(omegal)));
for i = 1:numel(omegal)
    if omegal(i) == omegar(i)
        vx(i) = omegal(i) * 0.2;
        vy(i) = 0;
        omegaz(i) = 0;
        continue
    end
    iomegal = omegal(i);
    iomegar = omegar(i);
    disp("sim for omegal = "+iomegal+", omegar = "+iomegar);
    out = sim("trkwhlmdldynf.slx");
    vx = out.vel.vx.Data(end);
    vy = out.vel.vy.Data(end);
    omegaz = out.vel.omegaz.Data(end);
    vx(i) = vx;
    vy(i) = vy;
    omegaz(i) = omegaz;
end

vsl = omegal.*0.2;
vsr = omegar.*0.2;
yc = vx./omegaz;
yl = (vx-vsl)./omegaz;
yr = (vx-vsr)./omegaz;
xc = -vy./omegaz;

abs((vsl-vsr).*(vsl+vsr));
abs((vsl-vsr)./(vsl+vsr));

figure("Position",[100,500,500,400]);
surf(omegal,omegar,vx);
xlabel("omegal");
ylabel("omegar");
title("vx");
exportgraphics(gca,"pic\dynf_vx.png");
figure("Position",[600,500,500,400]);
surf(omegal,omegar,vy);
xlabel("omegal");
ylabel("omegar");
title("vy");
exportgraphics(gca,"pic\dynf_vy.png");
figure("Position",[1100,500,500,400]);
surf(omegal,omegar,omegaz);
xlabel("omegal");
ylabel("omegar");
title("omegaz");
exportgraphics(gca,"pic\dynf_omegaz.png");

save("sspeverifysim.mat",...
    "omegallist","omegarlist","omegalgrid","omegargrid","vx","vy","omegaz");
