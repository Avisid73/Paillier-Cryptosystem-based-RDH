function plotHist(img_u8, ttl)
% plotHist  -- robust histogram display for RDHEI pipeline images.
%
% Why not just imhist()?
%   (1) imhist() on a double image assumes [0,1] range -> blank plot if
%       values are 0..255.
%   (2) Even on uint8, when one bin dwarfs all others by >20x (e.g., the
%       plain-stego histogram where half the non-seed pixels sit at 0),
%       imhist's default y-axis auto-scaling makes every other bin <1
%       pixel tall and the plot looks empty.
%
% This helper:
%   - computes the 256-bin count with imhist(uint8(.))
%   - draws with bar() (no display-scale trickery)
%   - auto-caps y at the 97th percentile of non-zero bins so dominant
%     bins don't hide the rest
%   - annotates any off-chart bin so the user knows it's not a bug

if ~isa(img_u8,'uint8'), img_u8 = uint8(img_u8); end

h = imhist(img_u8);                % 256x1 counts
x = 0:255;

% determine y cap from distribution of non-zero bin heights
nz  = h(h > 0);
if isempty(nz)
    ycap = 1;
else
    % 97th percentile -> ignore single huge spike
    ycap = max(10, double(prctile(nz, 97)) * 1.1);
end

bar(x, h, 'FaceColor', [0.12 0.47 0.71], 'EdgeColor', 'none', 'BarWidth', 1);
xlim([0 255]);
ylim([0 ycap]);
title(ttl);
xlabel('pixel value'); ylabel('count');

% annotate any bars clipped off the top
clipped = find(h > ycap);
if ~isempty(clipped)
    for k = 1:min(3, numel(clipped))
        b = clipped(k);
        text(double(b), ycap*0.95, sprintf('bin %d = %d', b, h(b)), ...
            'Color','r','FontSize',8,'HorizontalAlignment','left', ...
            'VerticalAlignment','top','Rotation',0);
    end
end
end
