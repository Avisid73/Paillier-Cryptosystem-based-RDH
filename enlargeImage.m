function [SI, seedMask] = enlargeImage(I)
% Scale image from P x Q to (2P-1) x (2Q-1) by inserting one row / column
% of zeros between every pair of original rows / columns (Section 4 of
%   Seed pixels (original pixel positions) are marked.
%
%  Input  : I (P x Q, uint8)
%  Output : SI       -- double image of size (2P-1) x (2Q-1)
%           seedMask -- logical, true at seed pixel positions

I = double(I);
[P,Q] = size(I);
SI = zeros(2*P-1, 2*Q-1);
SI(1:2:end, 1:2:end) = I;          % seed (original) pixels

% The paper fills non-seed pixels at the start as 0; the central non-seed
% of each 3x3 block later gets the integer average of the 4 corner seeds.
seedMask = false(size(SI));
seedMask(1:2:end, 1:2:end) = true;

% Fill centre of every 3x3 block with floor(avg of 4 corner seeds).
for r = 2:2:size(SI,1)-1
    for c = 2:2:size(SI,2)-1
        corners = [SI(r-1,c-1) SI(r-1,c+1) SI(r+1,c-1) SI(r+1,c+1)];
        SI(r,c) = floor(mean(corners));
    end
end
end
