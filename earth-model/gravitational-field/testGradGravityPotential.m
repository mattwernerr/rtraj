function testGradGravityPotential()
%TESTGRADGRAVITYPOTENTIAL Check GradGravityPotential against finite differences.
%   The output of GradGravityPotential is grad(V) in ECF spherical
%   coordinates, where V is the spherical harmonic potential
%
%       V = (GM/r) * (1 + sum_{n>=2} sum_{m<=min(n,M)} (Req/r)^n Pbar_nm(cos th) (C_nm cos(m lam) + S_nm sin(m lam))).
%
%   Each component of the gradient is compared with a central finite
%   difference of V. Run with MATLAB or Octave:
%
%       testGradGravityPotential()

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'math'), fullfile(here, '..', '..', 'check'));

GM = 3.986004418e14;   % Earth gravitational parameter [m3/s2]
Req = 6378137.0;       % Equatorial radius [m]
N = 4;                 % Maximum degree
M = 4;                 % Maximum order

% Coefficient ordering expected by GradGravityPotential: degree-major,
% with orders m = 0..min(n, M) for each degree n = 2..N.
n = zeros(0, 1);
m = zeros(0, 1);
for nn = 2:N
    for mm = 0:min(nn, M)
        n(end+1, 1) = nn; %#ok<AGROW>
        m(end+1, 1) = mm; %#ok<AGROW>
    end
end
L = numel(n);

% Deterministic, non-trivial coefficients (order 1e-3, the size of J2).
C = 1e-3 * sin(1:L)';
S = 1e-3 * cos(2.3 * (1:L))';

radii = [1.1, 1.6] * Req;
colats = [0.4, 1.2, 2.1];
lons = [0.3, 2.5, 5.9];

hr = 1;        % Radial step [m]
ht = 1e-6;     % Colatitude step [rad]
hl = 1e-6;     % Longitude step [rad]

for r = radii
    for th = colats
        for lam = lons
            g = GradGravityPotential(GM, Req, n, m, C, S, r, th, lam);

            dVdr = (potentialValue(GM, Req, n, m, C, S, r + hr, th, lam) ...
                - potentialValue(GM, Req, n, m, C, S, r - hr, th, lam)) / (2 * hr);
            dVdth = (potentialValue(GM, Req, n, m, C, S, r, th + ht, lam) ...
                - potentialValue(GM, Req, n, m, C, S, r, th - ht, lam)) / (2 * ht);
            dVdlam = (potentialValue(GM, Req, n, m, C, S, r, th, lam + hl) ...
                - potentialValue(GM, Req, n, m, C, S, r, th, lam - hl)) / (2 * hl);

            expected = [dVdr; dVdth / r; dVdlam / (r * sin(th))];

            % Tolerance: 1e-7 of the Newtonian gravity GM/r^2.
            tol = 1e-7 * GM / r^2;
            assert(all(abs(g - expected) <= tol), ...
                'GradGravityPotential differs from finite differences at r = %g, colatitude = %g, longitude = %g', ...
                r, th, lam);
        end
    end
end

disp('testGradGravityPotential: all finite-difference checks passed')
end

function V = potentialValue(GM, Req, n, m, C, S, r, th, lam)
%POTENTIALVALUE Evaluate the truncated spherical harmonic potential V.
V = 1;
for k = 1:numel(n)
    Pnm = computefnLegendre(n(k), cos(th), m(k));
    V = V + (Req / r)^n(k) * Pnm * (C(k) * cos(m(k) * lam) + S(k) * sin(m(k) * lam));
end
V = GM / r * V;
end
