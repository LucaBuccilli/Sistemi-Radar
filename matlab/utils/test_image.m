function img = test_image(n)
%TEST_IMAGE Immagine di prova sintetica in [0,1]: sfumatura, rettangolo, disco e righe oblique.
%   Si usa un'immagine sintetica per non dipendere da file esterni o da toolbox.
if nargin < 1, n = 256; end
[x, y] = meshgrid((0:n-1)/n, (0:n-1)/n);
img = 0.3 + 0.3*x;                                        % sfumatura (basse frequenze)
img(x > 0.15 & x < 0.45 & y > 0.2 & y < 0.5) = 0.9;       % rettangolo (bordi netti)
img((x-0.7).^2 + (y-0.35).^2 < 0.02) = 0.1;               % disco
stripes = 0.5 + 0.4*sin(2*pi*24*(x + y));                 % tessitura a righe oblique
mask = y > 0.65;
img(mask) = stripes(mask);
end
