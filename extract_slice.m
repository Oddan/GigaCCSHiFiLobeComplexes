load distal_fan_plane_N2S_processed

slice_dir = 'y'; % 'x' or 'y'


% create map from G actual cells to logical cells
remap = zeros(prod(G.cartDims),1);
remap(G.cells.indexMap) = 1:length(G.cells.indexMap);

% determine logical indices of cells in slice
lgrid = zeros(G.cartDims);

if slice_dir == 'y'
    slice_ix = 60;
    lgrid(:, slice_ix, :) = 1;
else
    slice_ix = 100;
    lgrid(slice_ix, :, :) = 1;
end

lgrid = logical(lgrid);

% determine actual indices of cells to keep in G
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

figure
plotCellData(Gs, log10(rocks.perm/1e3), 'edgecolor', 'none'); colorbar
if slice_dir == 'y'
    view(0,0);
else
    view(90,0);
end

set(gcf, 'Position', [471, 315, 2151, 797]);

% identify boundaries between permeable and practically impermeable layers
% by looking at a vertical column in the middle of the slice
lgrid = lgrid * 0;
remap_slice = zeros(prod(Gs.cartDims),1);
remap_slice(Gs.cells.indexMap) = 1:length(Gs.cells.indexMap);

if slice_dir == 'y'
    col_ix = 150;
    lgrid(col_ix, slice_ix, :) = 1;
    
else
    col_ix = round(Gs.cartDims(2)/2);
    lgrid(slice_ix, col_ix, :) = 1;
end
lgrid = logical(lgrid);
col_cells = remap_slice(lgrid(:));
col_cells = col_cells(col_cells > 0);

% find the two internal z-layers with very low permeability
thres = -4; % in log10(mD)
low_perm_layers = find(log10(rocks.permz(col_cells)/1e3) < thres);

% extract zone between low perm layers as separate grid
%keep_zone = (low_perm_layers(1)+4 : low_perm_layers(2)-1)
%keep_zone = (20:25);

lgrid = lgrid * 0;

if slice_dir == 'y'
    slice_ix = 60;
    lgrid(:, slice_ix, :) = 1;
else
    slice_ix = 100;
    lgrid(slice_ix, :, :) = 1;
end
lgrid(:, :, 1:(low_perm_layers(1))) = 0;
lgrid(:,:, low_perm_layers(2):end) = 0;

lgrid = logical(lgrid);
midlayer_ixs = remap(lgrid(:));
midlayer_ixs = midlayer_ixs(midlayer_ixs > 0);

Gmid = extractSubgrid(G, midlayer_ixs);
Gmid = computeGeometry(Gmid);
rocks_mid = rock;
rocks_mid.poro = rock.poro(midlayer_ixs);
rocks_mid.perm = rock.perm(midlayer_ixs, :);
rocks_mid.ntg = rock.ntg(midlayer_ixs);
rocks_mid.facies = rock.facies(midlayer_ixs);
rocks_mid.permz = rock.permz(midlayer_ixs);

figure
plotCellData(Gmid, log10(rocks_mid.perm/1e3), 'edgecolor', 'black'); colorbar
if slice_dir == 'y'
    view(0,0);
else
    view(90,0);
end

save('grids', 'Gs', 'Gmid', 'rocks', 'rocks_mid');
