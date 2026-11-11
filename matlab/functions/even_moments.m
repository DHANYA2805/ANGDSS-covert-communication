function S = even_moments(S, x, orders)
%EVEN_MOMENTS  Accumulate power sums of real samples for moment estimation.
%   S = even_moments(S, x, orders) adds sum(x), sum(x.^k) for the even
%   orders given (e.g. 2:2:20) to the running struct S. Start with S = [].
%   Use standardized_moments(S) at the end. Works chunk by chunk, so very
%   long signals (1e6 bits x N = 256) never have to be held in memory.
x = x(:);
if isempty(S)
    S.n = 0; S.s1 = 0; S.orders = orders; S.sk = zeros(size(orders));
end
S.n  = S.n + numel(x);
S.s1 = S.s1 + sum(x);
x2 = x.^2;
if isequal(S.orders, [2 20])           % fast path for the 20th-moment detector
    x4 = x2.^2; x16 = (x4.^2).^2;
    S.sk = S.sk + [sum(x2), sum(x16 .* x4)];
    return;
end
p  = x2;
o  = 2;
for k = 1:max(orders)/2
    j = find(S.orders == o, 1);
    if ~isempty(j), S.sk(j) = S.sk(j) + sum(p); end
    p = p .* x2;
    o = o + 2;
end
end
