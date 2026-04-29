function [recovered, W] = extractData(decStego, seedMask, P, Q, sideInfo) %#ok<INUSL,INUSD>
% EXTRACTDATA  Algorithm 2 
%
% Self-contained: reads the payload length from the first Llen central-
% non-seed LSBs, derives how many seeds were actually embedded, and walks
% blocks in the SAME FORWARD order used by embedData.  Stops after
% recovering exactly the embedded seed count so unused seeds are left
% untouched (critical for bit-exact image recovery when payload < capacity).

[H, Wd] = size(decStego);
img = double(decStego);

% same Llen formula as embed
Llen = max(16, ceil(log2(max(9*P*Q, 2))) + 1);

% --- centres -----------------------------------------------------------
centres = zeros(0,2);
for br = 2:2:H-1
    for bc = 2:2:Wd-1
        centres(end+1,:) = [br bc]; %#ok<AGROW>
    end
end
nBlocks = size(centres,1);

% --- read length header -------------------------------------------------
centreBits = zeros(nBlocks,1);
for i = 1:nBlocks
    centreBits(i) = mod(img(centres(i,1), centres(i,2)), 2);
end
payloadLen   = bi2de(centreBits(1:Llen).', 'left-msb');
nCentreFree  = nBlocks - Llen;
nCentrePay   = min(nCentreFree, payloadLen);
nSeedBits    = payloadLen - nCentrePay;
padBits      = mod(nSeedBits, 2);
nU           = (nSeedBits + padBits) / 2;

W_centres = centreBits(Llen+1 : Llen+nCentrePay);

% --- block walk, stop after nU seed recoveries --------------------------
seedRecovered = false(size(img));
Uvals = zeros(nU,1);
uCount = 0;
done = false;
for br = 2:2:H-1
    if done, break; end
    for bc = 2:2:Wd-1
        if done, break; end
        seeds = [br-1 bc-1;
                 br-1 bc+1;
                 br+1 bc+1;
                 br+1 bc-1];
        nexts = [br-1 bc;
                 br   bc+1;
                 br+1 bc;
                 br   bc-1];

        for k = 1:4
            sr = seeds(k,1); sc = seeds(k,2);
            if seedRecovered(sr,sc), continue; end
            if uCount >= nU
                done = true; break;
            end

            nr = nexts(k,1); nc = nexts(k,2);
            byteVal = img(nr,nc);
            Qs  = floor(byteVal/10);
            dir = mod(byteVal,10);

            s = 256*Qs + img(sr,sc);
            U = mod(s,4);
            uCount = uCount + 1;
            Uvals(uCount) = U;

            if dir == 0
                X = (s - U)/2;                        % (3b)
            else
                X = (s - U - 512)/2 + 1;              % (3a)
            end
            img(sr,sc) = X;
            seedRecovered(sr,sc) = true;
        end
    end
end

% seed bits, forward embedding order
if nU > 0
    bitsFromSeeds = reshape(de2bi(Uvals,2,'left-msb').', [], 1);
else
    bitsFromSeeds = [];
end
if padBits > 0 && ~isempty(bitsFromSeeds)
    bitsFromSeeds = bitsFromSeeds(1:end-padBits);
end

W = [W_centres(:); bitsFromSeeds(:)];
if numel(W) >= payloadLen
    W = W(1:payloadLen);
end

% drop non-seeds: original cover was the (odd, odd) subgrid of SI
recoveredSeeds = img(1:2:end, 1:2:end);
recovered = uint8(recoveredSeeds);
end
