function M = paillierDecrypt(C, N, lambda, mu)
% PAILLIERDECRYPT  Paillier decryption per-pixel.
%
%   m = paillierDecrypt(c, N, lambda, mu)
%
% Plaintext:  m = L(c^lambda mod N^2) * mu   (mod N)
% where L(u) = (u-1)/N.
%
% Paper's keys: (N=7387, lambda=3608, mu=7324).

C = double(C);
N2 = N * N;
[rows, cols] = size(C);
M = zeros(rows, cols);

for i = 1:numel(C)
    u  = modpow(C(i), lambda, N2);     % u = c^lambda mod N^2
    Lv = (u - 1) / N;                  % integer by construction
    M(i) = mod(Lv * mu, N);
end
end

function y = modpow(base, exp, m)
y = 1;
base = mod(base, m);
while exp > 0
    if mod(exp,2) == 1
        y = mod(y * base, m);
    end
    exp = floor(exp/2);
    base = mod(base * base, m);
end
end
