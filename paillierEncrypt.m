function C = paillierEncrypt(M, N, g, r)
% PAILLIERENCRYPT  Paillier encryption per-pixel.
%
%   C = paillierEncrypt(M, N, g)        -- uses a fixed deterministic r
%   C = paillierEncrypt(M, N, g, r)     -- user-supplied r (scalar or matrix)
%
% Ciphertext:   c = g^m * r^N  (mod N^2)
% Plaintext m must be integer in [0, N-1].
% r must be a random integer in [1, N-1] with gcd(r,N)=1.
%
% This implements the paper's Paillier cryptosystem with the keys
% (N=7387, g=75).  Each pixel becomes an integer in [0, N^2-1] =
% [0, 54,567,768], which fits in a double (exact up to 2^53 ~= 9e15).
%


if nargin < 4
    r = [];
end

M  = double(M);
N2 = N * N;
assert(all(M(:) >= 0 & M(:) < N), ...
       'paillierEncrypt: plaintext out of [0,N-1]=[0,%d]', N-1);

[rows, cols] = size(M);
C = zeros(rows, cols);

% Deterministic r from a PRNG (keyed elsewhere).  
if isempty(r)
    s = RandStream('mt19937ar','Seed',uint32(987654321));
    r = randi(s, [2 N-1], rows, cols);
    % ensure gcd(r,N)=1.  With N=7387=7*1055? (actually 7387=83*89) this
    % is almost always true; fix if not.
    for i = 1:numel(r)
        while gcd(r(i), N) ~= 1
            r(i) = mod(r(i)+1, N-1) + 2;
        end
    end
end

% c = (g^m * r^N) mod N^2
for i = 1:numel(M)
    gm = modpow(g, M(i), N2);
    rn = modpow(r(i), N, N2);
    C(i) = mod(gm * rn, N2);
end
end

% -------- helper: modular exponentiation by squaring --------
function y = modpow(base, exp, m)
% y = base^exp mod m, using double-precision arithmetic.
% Safe while m < 2^26 (so m*m < 2^52 stays exact in double).
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
