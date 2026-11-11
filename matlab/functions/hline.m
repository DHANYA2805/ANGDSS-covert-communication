function h = hline(y, style, label)
%HLINE  Horizontal reference line (portable replacement for yline).
xl = xlim;
h  = plot(xl, [y y], style, 'LineWidth', 1.2);
if nargin > 2 && ~isempty(label)
    text(xl(1) + 0.02 * diff(xl), y, ['  ' label], 'VerticalAlignment', 'bottom', ...
        'Color', get(h, 'Color'), 'FontSize', 9);
end
end
