function ax = mod_prog(epsilon, transfer, varargin)
if nargin == 2
    figure;
    ax = newplot;
    epsilon(transfer == 0) = [];
    transfer(transfer == 0) = [];
    semilogy(ax, epsilon, abs(transfer), 'o');
    ax.Title.String = 'Line Search for Sign Change';
    ax.XLabel.String = 'Epsilon Value';
    ax.YLabel.String = 'Unresolved Transfer Term Magnitude';
else
    ax = varargin{1};
    epsilon(transfer == 0) = [];
    transfer(transfer == 0) = [];
    semilogy(ax, epsilon, abs(transfer), 'o');
end
drawnow;