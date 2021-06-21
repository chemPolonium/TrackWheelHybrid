omegallist = 19:0.2:21;
omegarlist = 19:0.2:21;
[omegalgrid,omegargrid] = meshgrid(omegallist,omegarlist);

% [vx,vy,omegaz] = arrayfun(@finalState,omegalgrid,omegargrid);
[vx,vy,omegaz] = deal(zeros(size(omegalgrid)));
for i = 1:numel(omegalgrid)
    disp(i);
%     if omegalgrid(i) == omegargrid(i)
%         vx(i) = omegalgrid(i) * 0.3;
%         vy(i) = 0;
%         omegaz(i) = 0;
%         continue
%     end
    omegal = omegalgrid(i);
    omegar = omegargrid(i);
    out = sim("trackwheelmodel.slx");
    vx(i) = out.vx.Data(end);
    vy(i) = out.vy.Data(end);
    omegaz(i) = out.omegaz.Data(end);
end

surf(omegalgrid,omegargrid,vx);
title("vx");
exportgraphics(gca,"pic\vx.png");
surf(omegalgrid,omegargrid,vy);
title("vy");
exportgraphics(gca,"pic\vy.png");
surf(omegalgrid,omegargrid,omegaz);
title("omegaz");
exportgraphics(gca,"pic\omegaz.png");

save("sspeverifysim.mat",...
    "omegallist","omegarlist","omegalgrid","omegargrid","vx","vy","omegaz");

% function [vx,vy,omegaz] = finalState(omegal,omegar)
%     assignin("base","omegal",omegal);
%     assignin("base","omegar",omegar);
%     out = sim("trackwheelmodel.slx");
%     vx = out.vx.Data(end);
%     vy = out.vy.Data(end);
%     omegaz = out.omegaz.Data(end);
% end