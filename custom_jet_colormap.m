function cmap = custom_jet_colormap(n)
% CUSTOM_JET_COLORMAP  Mimics the red→orange→yellow→green→cyan→blue colormap
%   cmap = custom_jet_colormap(n)  returns an n-by-3 colormap matrix.
%   Default n = 256 if not specified.

% The specific purpose of this colormap is to represent simulation results using
% the same (or similar) color scheme as the default in ResInsight.

if nargin < 1
    n = 256;
end

% Key color stops: [R, G, B] from top (high) to bottom (low)
% Reversed so index 1 = bottom (deep blue), index end = top (red)
% stops = [
%     0.60,  0.00,  0.80;   % deep blue/indigo  (bottom)
%     0.00,  0.20,  1.00;   % blue
%     0.00,  0.80,  1.00;   % cyan
%     0.00,  1.00,  0.20;   % green
%     1.00,  1.00,  0.00;   % yellow
%     1.00,  0.50,  0.00;   % orange
%     1.00,  0.00,  0.00;   % red               (top)
% ];

stops = [
    0.00,  0.00,  0.55;   % dark blue          (bottom)
    0.00,  0.20,  1.00;   % blue
    0.00,  0.80,  1.00;   % cyan
    0.00,  1.00,  0.20;   % green
    1.00,  1.00,  0.00;   % yellow
    1.00,  0.50,  0.00;   % orange
    1.00,  0.00,  0.00;   % red                (top)
];

% Interpolate across n levels
x_stops = linspace(0, 1, size(stops, 1));
x_query = linspace(0, 1, n);

cmap = zeros(n, 3);
for ch = 1:3
    cmap(:, ch) = interp1(x_stops, stops(:, ch), x_query, 'pchip');
end

% Clamp to [0, 1]
cmap = max(0, min(1, cmap));

end
