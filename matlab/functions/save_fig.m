function save_fig(fh, cfg, name)
%SAVE_FIG  Save a figure as PNG (and .fig in MATLAB) into cfg.outdir.
if ~cfg.saveFigs, return; end
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end
set(fh, 'PaperPositionMode', 'auto');
print(fh, fullfile(cfg.outdir, [name '.png']), '-dpng', '-r150');
if ~exist('OCTAVE_VERSION', 'builtin')
    try, savefig(fh, fullfile(cfg.outdir, [name '.fig'])); catch, end
end
end
