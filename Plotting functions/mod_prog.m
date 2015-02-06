function ax = mod_prog(epsilon, transfer, varargin)
if nargin == 2
    ax = newplot;
else
    ax = varargin{1};
end
plot(ax, epsilon, transfer);
ax.Title.String = 'Line Search for Sign Change';
ax.XLabel.String = 'Epsilon Value';
ax.YLabel.String = 'Unresolved Transfer Term';
drawnow;