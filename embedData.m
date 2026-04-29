function [stego, sideInfo] = embedData(SI, seedMask, W, P, Q) %#ok<INUSL>
% EMBEDDATA  Algorithm 1 of Singh & Singh (Sec. 4.1), one-embedding-per-

%
% Inputs :
%   SI        -- enlarged cover image (2P-1) x (2Q-1), double in [0,255]
%   seedMask  -- logical mask of seed pixel positions (kept for API
%                symmetry; not needed internally)
%   W         -- payload bit vector (0/1)
%   P,Q       -- dimensions of the ORIGINAL cover image
%
% Outputs :
%   stego     -- plaintext stego image, strictly in [0,255]  (double)
%   sideInfo  -- struct echoing quantities useful for diagnostics;
%                the extractor does NOT require sideInfo because it reads
%                the payload length from the length-header bits.
%
% Core math (from the paper) --------------------------------------------
%   (2a)  seed odd  : s = 2*(X-1) + U + 512
%   (2b)  seed even : s = 2*X     + U                U in {0,1,2,3}
%
% FOB split --------------------------------------------------------------
%   After (2a)/(2b), s can reach ~1023.  We split s = 256*Qs + Rs.  Rs
%   stays in the seed; Qs is decimal-packed into the "next non-seed":
%         non_seed_byte = 10*Qs + dir_code,    dir_code in {0,1,2}
%   dir_code = 0  -> (2b) was used
%   dir_code = 1  -> (2a) + next non-seed is to the right of, or above, the seed
%   dir_code = 2  -> (2a) + next non-seed is to the left of, or below, the seed
%
% Length header ----------------------------------------------------------
%   Paper uses Llen = log2(P)+log2(Q), which is too small for the scheme's
%   capacity (e.g. P=Q=64 gives Llen=12 -> max 4096 representable, but
%   capacity is ~11000 bits).  We use
%       Llen = max(16, ceil(log2(9*P*Q))+1)
%   which is comfortably enough to index the true capacity.  This is a
%   corrective deviation from the paper.

[H, Wd] = size(SI);
stego   = SI;                            % work in double throughout

% --- length header size (corrected) -------------------------------------
Llen  = max(16, ceil(log2(max(9*P*Q, 2))) + 1);
Lbits = de2bi(numel(W), Llen, 'left-msb');

% --- central non-seed positions (block centres) ------------------------
centres = zeros(0,2);
for br = 2:2:H-1
    for bc = 2:2:Wd-1
        centres(end+1,:) = [br bc]; %#ok<AGROW>
    end
end
nBlocks = size(centres,1);
assert(numel(Lbits) <= nBlocks, ...
    'Image too small to hold length header: need %d centres, have %d', ...
    numel(Lbits), nBlocks);

% --- Step 1: embed length header via even-parity LSB on first Llen centres
% Boundary-safe: if adding 1 would push a 255-valued pixel to 256, we
% subtract 1 instead (both toggle the LSB and both are reversible since
% the extractor only reads LSB).  This is the "reserved boundary bin 255"
% safeguard; bin 0 never needs protection because LSB(0)=0 and if target
% is 1 we add 1 -> 1 which is in range.
for i = 1:Llen
    r = centres(i,1); c = centres(i,2);
    if mod(stego(r,c),2) ~= Lbits(i)
        if stego(r,c) >= 255
            stego(r,c) = stego(r,c) - 1;
        else
            stego(r,c) = stego(r,c) + 1;
        end
    end
end

% --- Step 2: embed leading payload bits in the remaining centres --------
nCentreFree = nBlocks - Llen;
nCentrePay  = min(nCentreFree, numel(W));
for i = 1:nCentrePay
    r = centres(Llen+i,1); c = centres(Llen+i,2);
    b = W(i);
    if mod(stego(r,c),2) ~= b
        if stego(r,c) >= 255
            stego(r,c) = stego(r,c) - 1;
        else
            stego(r,c) = stego(r,c) + 1;
        end
    end
end
remaining = W(nCentrePay+1:end);
remaining = remaining(:);

% --- pair the remaining bits into base-4 U values -----------------------
padBits = 0;
if mod(numel(remaining),2) == 1
    remaining(end+1) = 0;
    padBits = 1;
end
if isempty(remaining)
    Uvals = [];
else
    Uvals = bi2de(reshape(remaining,2,[]).', 'left-msb');
end

% --- Step 3: block walk, one embedding per seed -------------------------
processedSeed = false(size(stego));
u_idx = 1;
for br = 2:2:H-1
    for bc = 2:2:Wd-1
        seeds = [br-1 bc-1;          % TL
                 br-1 bc+1;          % TR
                 br+1 bc+1;          % BR
                 br+1 bc-1];         % BL
        nexts = [br-1 bc;            % TL -> top-centre
                 br   bc+1;          % TR -> right-centre
                 br+1 bc;            % BR -> bottom-centre
                 br   bc-1];         % BL -> left-centre

        for k = 1:4
            sr = seeds(k,1); sc = seeds(k,2);
            if processedSeed(sr,sc),        continue; end
            if u_idx > numel(Uvals),        break;    end
            U = Uvals(u_idx); u_idx = u_idx + 1;

            X = stego(sr,sc);                     % double, in [0,255]
            if mod(X,2) == 1
                s = (X-1)*2 + U + 512;             % (2a)
                used2a = true;
            else
                s = X*2 + U;                       % (2b)
                used2a = false;
            end

            % FOB split
            Qs = floor(s/256);
            Rs = mod(s,256);
            stego(sr,sc) = Rs;

            % direction code on the next non-seed
            nr = nexts(k,1); nc = nexts(k,2);
            dir = 0;
            if used2a
                drow = nr - sr;   % +1 below, -1 above
                dcol = nc - sc;   % +1 right, -1 left
                if drow > 0 || dcol < 0          % below or left
                    dir = 2;
                elseif drow < 0 || dcol > 0      % above or right
                    dir = 1;
                end
            end
            stego(nr,nc) = 10*Qs + dir;            % replaces prior value

            processedSeed(sr,sc) = true;
        end
        if u_idx > numel(Uvals), break; end
    end
    if u_idx > numel(Uvals), break; end
end

% Invariant: no pixel out of range (no clipping used, all by construction).
assert(all(stego(:) >= 0 & stego(:) <= 255), ...
    'FOB violation: stego pixel out of [0,255] (this should be impossible).');

sideInfo = struct('padBits', padBits, ...
                  'nPayload', numel(W), ...
                  'nCentrePay', nCentrePay, ...
                  'Llen', Llen);
end
