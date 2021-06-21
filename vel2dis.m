function [x,y] = vel2dis(u,v,r,t,cx,cy)

psi = r.*t;
x = u.*sinc(psi).*t - v.*cosc(psi).*t + cx.*cos(psi) - cy.*sin(psi);
y = v.*sinc(psi).*t + u.*cosc(psi).*t + cx.*sin(psi) + cy.*cos(psi);

end

function y = sinc(x)
y = ones(size(x));
logi = abs(x) > 0.01;
y(logi) = sin(x(logi))./x(logi);
end

function y = cosc(x)
y = zeros(size(x))+x;
logi = abs(x) > 0.01;
y(logi) = (1-cos(x(logi)))./x(logi);
end