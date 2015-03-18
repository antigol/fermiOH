function y = approx(ef, r)
% E_F = ef hbar omega
% r[m] = r / sqrt(nu) ; nu = mass omega / (2 hbar)
% y = density in units of nu^(3/2)

y = 16 / (3 * pi) * (ef - r.^2).^(3/2)  .* (r <= sqrt(ef));
