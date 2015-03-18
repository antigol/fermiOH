# http://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#N-dimensional_harmonic_oscillator
function result = fermiOH(nmax, r)
% nu = mass omega / (2 hbar)   [m^-2]

% n = 2k+l ; E = hbar omega (n + 3/2) (si on ignore le 3/2, E = hbar omega n)

%r = linspace(0, 10, 1000); % in units of nu^(-1/2)
r = reshape(r, 1, numel(r));
result = zeros(size(r));   % in units of nu^(3/2)

for k=0:nmax/2,
for l=0:nmax-2*k,
	P = LaguerreGen(k, l + 0.5);
	
	f = prod([(k:-1:1) ones(1, l+1)] ./ [2*k+2*l+1:-2:1]);
	result = result + (2*l+1) * sqrt(2/pi) * 2^(k+2*l+3) * f * r.^(2*l) .* exp(-2*r.^2) .* polyval(P, 2*r.^2).^2;
end
end
