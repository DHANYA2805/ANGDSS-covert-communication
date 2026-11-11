function z = despread(r, pn)
%DESPREAD  Correlate received chips with the binary PN code.
%   z = despread(r, pn) returns one decision variable per symbol (column).
%   The receiver uses only the +/-1 code; it never needs the mask |U|,|V|.
N    = size(pn, 1);
nSym = numel(r) / N;
R    = reshape(r, N, nSym);
if size(pn, 2) == 1
    z = (pn.' * R).';
else
    z = sum(R .* pn, 1).';
end
end
