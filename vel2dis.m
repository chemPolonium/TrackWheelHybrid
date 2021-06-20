function [x,y] = vel2dis(u,v,r,t,cx,cy)

psi = r.*t;
x = u.*sinc(psi).*t + v.*cosc(psi).*t + cx.*cos(psi) - cy.*sin(psi);
y = v.*sinc(psi).*t - u.*cosc(psi).*t + cx.*sin(psi) + cy.*cos(psi);

end

function y = sinc(x)
y = ones(size(x));
y(x~=0) = sin(x(x~=0))./x(x~=0);
end

function y = cosc(x)
y = zeros(size(x));
y(x~=0) = (cos(x(x~=0))-1)./x(x~=0);
end