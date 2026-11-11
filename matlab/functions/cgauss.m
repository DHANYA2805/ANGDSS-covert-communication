function w = cgauss(n, v)
%CGAUSS  Circularly-symmetric complex Gaussian noise, total variance v.
%   w = cgauss(n, v): n x 1 samples, E|w|^2 = v. v may be a scalar or n x 1.
w = sqrt(v / 2) .* (randn(n, 1) + 1j * randn(n, 1));
end
