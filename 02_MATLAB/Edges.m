
negtv = CDASIMG;
negtv(negtv < 0) = 0;
negtv = abs(negtv);

%compare edge detection algorithms
algo_type = ["Sobel" "Prewitt" "Roberts"];
tiledlayout(2, length(algo_type));

for ii = algo_type
nexttile;
img = edge(CDASIMG,ii);
imagesc(img);
title(ii)
subtitle("Positive and Negative")
end

for ii = algo_type
nexttile;
img = edge(negtv,ii);
imagesc(img);
title(ii)
subtitle("Positive only")
colormap(pucolors.viridis)
end