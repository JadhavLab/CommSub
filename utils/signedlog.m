function y = signedlog(x)

nanlocs = isnan(x);
x(nanlocs) = 0;
s = sign(x);
x = abs(x);
x = log10(x + 10) - 1;
y = x .* s;
y(nanlocs) = nan;
