basefunc = @(x,m)1-(1-exp(-m.*x))./(m.*x);

fplot(@(x)basefunc(x,1:5),[0 20])
