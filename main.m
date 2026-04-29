%% main.m -- 
%
%  Uses the Paillier cryptosystem with the paper's exact keys:
%     Encryption:  (N=7387, g=75)
%     Decryption:  (lambda=3608, mu=7324)
%
%  Pipeline:
%    cover -> enlarge -> plaintext embed (+FOB split to [0,255])
%          -> Paillier-encrypt each pixel (mod N^2 = 54,567,769)
%          -> transmit
%          -> Paillier-decrypt each pixel back to [0,255]
%          -> extract + reverse FOB -> original cover + 48 embedded bits

clear all; close all; clc;   %#ok<CLALL>
rng(2026);

fprintf('\n============================================================\n');
fprintf('  RDHEI with PAILLIER (keys: N=7387, g=75)\n');
fprintf('  Build 2026-04-18, binary payload, Baboon cover\n');
fprintf('============================================================\n\n');

% ---- 1. Load cover -------------------------------------------------
I = imread('/Users/sudhirsingh/Desktop/files/Baboon.tiff');
if size(I,3)==3, I = rgb2gray(I); end
I = uint8(I);
[P,Q] = size(I);
fprintf('Cover image size : %d x %d\n', P, Q);

% ---- 2. Payload: explicit 48-bit binary string ---------------------
secretBitString = '110101110001101101011001001111000101010011101101';
W = double(secretBitString - '0'); W = W(:);
fprintf('Secret payload   : %d bits\n', numel(W));
fprintf('Secret bits      : %s\n', secretBitString);

% ---- 3. Paillier keys (paper's toy keys, insecure but correct) ----
N      = 7387;
g      = 75;
lambda = 3608;
mu     = 7324;
N2     = N * N;
fprintf('Paillier keys    : N=%d, g=%d, lambda=%d, mu=%d\n', N, g, lambda, mu);
fprintf('                   N^2 = %d (ciphertext range [0, N^2-1])\n', N2);

% ---- 4. Plaintext embed + FOB --------------------------------------
[SI, seedMask]    = enlargeImage(I);
[stego, sideInfo] = embedData(SI, seedMask, W, P, Q);   % pixels in [0,255]
fprintf('Plaintext stego range (post-FOB): [%d, %d]\n', ...
        round(min(stego(:))), round(max(stego(:))));

% ---- 5. Paillier encrypt each pixel --------------------------------
% Paillier is probabilistic: c = g^m * r^N mod N^2 with random r.
% For a reproducible demo we use a fixed PRNG to derive r per pixel.
% Each pixel becomes an integer in [0, N^2-1] (won't fit in 8 bits).
fprintf('Encrypting with Paillier... ');
t0 = tic;
encStego = paillierEncrypt(stego, N, g);    % returns double matrix of ints
fprintf('done in %.2f s\n', toc(t0));
fprintf('Ciphertext range : [%d, %d]\n', min(encStego(:)), max(encStego(:)));

% ---- 6. Paillier decrypt -------------------------------------------
fprintf('Decrypting with Paillier... ');
t0 = tic;
decStego = paillierDecrypt(encStego, N, lambda, mu);
fprintf('done in %.2f s\n', toc(t0));
fprintf('Decrypted matches plaintext stego: %d\n', isequal(decStego, stego));

% ---- 7. Extract + recover ------------------------------------------
[recovered, Wout] = extractData(decStego, seedMask, P, Q, sideInfo);

% ---- 8. Verify ------------------------------------------------------
okImg  = isequal(recovered, I);
nBits  = min(numel(W), numel(Wout));
bitErr = sum(W(1:nBits) ~= Wout(1:nBits));
extractedStr = char('0' + Wout(1:numel(W)).');

fprintf('\n======= Verification =======\n');
fprintf('Image recovered perfectly : %d\n', okImg);
fprintf('Bit errors in payload     : %d / %d\n', bitErr, nBits);
fprintf('Original  bits : %s\n', secretBitString);
fprintf('Extracted bits : %s\n', extractedStr);
fprintf('Exact match    : %d\n', strcmp(secretBitString, extractedStr));

% ---- 9. Visualization -----------------------------------------------
% For the encrypted histogram we must map N^2-range ciphertexts down to
% uint8 for display.  We use mod 256, which is what the paper (implicitly)
% shows when drawing an "encrypted image".
I_u         = uint8(I);
stego_u     = uint8(stego);
encStego_u  = uint8(mod(encStego, 256));     % display-only reduction
recov_u     = uint8(recovered);

figure('Name','RDHEI / real Paillier','Color','w','Position',[50 50 1400 720]);
subplot(2,4,1); imshow(I_u);         title('Original cover');
subplot(2,4,2); imshow(stego_u);     title('Plain stego (post-FOB)');
subplot(2,4,3); imshow(encStego_u);  title('Paillier ciphertext');
subplot(2,4,4); imshow(recov_u);     title('Recovered cover');
subplot(2,4,5); plotHist(I_u,        'Hist: original');
subplot(2,4,6); plotHist(stego_u,    'Hist: plain marked');
subplot(2,4,7); plotHist(encStego_u, 'Hist: Paillier');
subplot(2,4,8); plotHist(recov_u,    'Hist: recovered');
saveas(gcf, 'rdhei_result.png');

% ---- 10. Chi-square uniformity -------------------------------------
h_enc    = imhist(encStego_u);
expected = numel(encStego_u) / 256;
chi2     = sum((h_enc - expected).^2) / expected;
fprintf('Paillier-mod-256 chi-square vs uniform = %.2f  (df=255, critical ~293)\n', chi2);
fprintf('Saved figure                           : rdhei_result.png\n');
