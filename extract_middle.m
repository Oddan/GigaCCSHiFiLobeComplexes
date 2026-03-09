mrstModule add mrst-gui

load distal_fan_plane_N2S_processed

% we will extract logical layers 115-227

% create map from G actual cells to logical cells
remap = zeros(prod(G.cartDims),1);
remap(G.cells.indexMap) = 1:length(G.cells.indexMap);

% determine actual indices of cells to keep in G
lgrid = zeros(G.cartDims);
lgrid(:, :, 115:227) = 1;
lgrid = logical(lgrid);

keep = remap(lgrid(:));
keep = keep(keep > 0);
Gs = extractSubgrid(G, keep);
Gs = computeGeometry(Gs);

% create the corresponding rock object
rocks = rock;
rocks.poro = rock.poro(keep);
rocks.perm = rock.perm(keep, :);
rocks.ntg = rock.ntg(keep);
rocks.facies = rock.facies(keep);
rocks.permz = rock.permz(keep);

% plot the resulting grid and rock properties
figure;
plotToolbar(Gs, rocks); camproj perspective;

% save grid
save('middle_grid_and_rock_compressed', 'Gs', 'rocks', '-v7.3');
