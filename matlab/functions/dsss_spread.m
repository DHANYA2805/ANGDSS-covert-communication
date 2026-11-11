function chips = dsss_spread(sym, pn)
%DSSS_SPREAD  Direct-sequence spreading of QPSK symbols.
%   chips = dsss_spread(sym, pn)
%   sym : nSym x 1 complex symbols
%   pn  : N x 1 (+/-1) short code repeated every symbol, or
%         N x nSym (+/-1) long code (a fresh N-chip code for every symbol)
%   chips : (N*nSym) x 1 column, symbol-major (N chips of symbol 1, then 2 ...)
sym  = sym(:).';
nSym = numel(sym);
N    = size(pn, 1);
if size(pn, 2) == 1
    C = pn * sym;                       % N x nSym (outer product)
else
    C = pn .* repmat(sym, N, 1);        % N x nSym
end
chips = C(:);
end
