function animateTop(t, R, opts)
% ANIMATETOP  3D animation of a spinning top whose tip rests on the floor.
%
%   animateTop(t, R)        plays the motion in a figure window.
%   animateTop(t, R, opts)  with options (all optional):
%       opts.video   = 'top.mp4'  also saves the animation as an MP4 file
%       opts.speed   = 1          playback speed (1/6 = six times slower)
%       opts.trail   = .5         length of the traced path of the stem (s)
%       opts.readout = {}         N x 1 cell array of text shown at each frame
%
%   t : N x 1 time vector (uniformly sampled, s)
%   R : 3 x 3 x N rotation matrices from the body frame to the space frame
%       (body z axis = symmetry axis of the top, tip at the origin)
%   The shape of the top is defined in topProfile.m.

if nargin < 3, opts = struct; end
if ~isfield(opts, 'video'), opts.video = ''; end
if ~isfield(opts, 'speed'), opts.speed = 1; end
if ~isfield(opts, 'trail'), opts.trail = .5; end
if ~isfield(opts, 'readout'), opts.readout = {}; end
N = numel(t); dt = t(2) - t(1);

% ---------------- Geometry (body frame) ----------------
P = topProfile();
nA = 60;
[Xb, Yb, Zb] = revolve(P.r, P.z, nA);
[Xs, Ys, Zs] = revolve([P.rs P.rs 0], [P.zs P.zs(end)], nA);
% Three red and three yellow sectors make the spin visible.
sector = mod(floor((0:nA)*6/nA), 2);
Cb = repmat(sector, numel(P.z), 1);
stemTip = [0; 0; P.zs(end)];
body = [Xb(:) Yb(:) Zb(:)].';               % for the shadow

% ---------------- Figure ----------------
fig = figure('Color', 'w', 'Position', [100 100 800 800], 'Name', 'Spinning top');
ax = axes(fig, 'Position', [0 0 1 1]); hold(ax, 'on'); axis(ax, 'equal', 'off');
lim = .09;
axis(ax, [-lim lim -lim lim 0 .08]);
view(ax, 35, 20); camtarget(ax, [0 0 .022]); camva(ax, 7.5);

% Floor with a 1 cm grid
[Fx, Fy] = meshgrid(linspace(-lim, lim, 2));
surf(ax, Fx, Fy, zeros(2), 'FaceColor', [.94 .95 .97], 'EdgeColor', 'none', 'FaceLighting', 'none');
for k = -lim:.01:lim
    plot3(ax, [k k], [-lim lim], [0 0] + 1e-5, 'Color', [.86 .87 .9]);
    plot3(ax, [-lim lim], [k k], [0 0] + 1e-5, 'Color', [.86 .87 .9]);
end

% Shadow: the body projected straight down onto the floor
shadow = surf(ax, Xb, Yb, zeros(size(Xb)) + 2e-5, 'FaceColor', [.2 .2 .25], 'FaceAlpha', .18, ...
    'EdgeColor', 'none', 'FaceLighting', 'none');

% Top (body and stem), drawn with hgtransform so only the matrix changes
tf = hgtransform(ax);
colormap(ax, [.80 .15 .15; .98 .80 .20]);
surf(ax, Xb, Yb, Zb, Cb, 'Parent', tf, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
    'AmbientStrength', .45, 'DiffuseStrength', .7, 'SpecularStrength', .5);
surf(ax, Xs, Ys, Zs, 'Parent', tf, 'FaceColor', [.25 .25 .3], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
caxis(ax, [0 1]);
trail = plot3(ax, nan, nan, nan, 'Color', [.15 .35 .75], 'LineWidth', 1.5);
light(ax, 'Position', [1 -1 2], 'Style', 'infinite'); light(ax, 'Position', [-1 1 .5], 'Style', 'infinite');
label = annotation(fig, 'textbox', [.04 .78 .7 .18], 'String', '', 'EdgeColor', 'none', ...
    'FontSize', 15, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'Color', [.15 .15 .2]);

% ---------------- Playback / recording ----------------
if ~isempty(opts.video)
    vw = VideoWriter(opts.video, 'MPEG-4'); vw.FrameRate = round(opts.speed/dt); vw.Quality = 95; open(vw);
end
nTrail = max(1, round(opts.trail/dt));
tipPath = zeros(3, N);
for k = 1:N
    Rk = R(:, :, k);
    M = eye(4); M(1:3, 1:3) = Rk; set(tf, 'Matrix', M);
    w = Rk*body;
    set(shadow, 'XData', reshape(w(1, :), size(Xb)), 'YData', reshape(w(2, :), size(Xb)));
    tipPath(:, k) = Rk*stemTip;
    i0 = max(1, k - nTrail);
    set(trail, 'XData', tipPath(1, i0:k), 'YData', tipPath(2, i0:k), 'ZData', tipPath(3, i0:k));
    if ~isempty(opts.readout), set(label, 'String', opts.readout{k}); end
    drawnow;
    if isempty(opts.video)
        pause(max(0, dt/opts.speed - .002));
        if ~isvalid(fig), return; end
    else
        writeVideo(vw, getframe(fig));
    end
end
if ~isempty(opts.video), close(vw); fprintf('Saved %s\n', opts.video); end
end

function [X, Y, Z] = revolve(r, z, n)
% Surface of revolution about the z axis.
a = linspace(0, 2*pi, n + 1);
X = r(:)*cos(a); Y = r(:)*sin(a); Z = repmat(z(:), 1, n + 1);
end
