function plotVerticalProfile(Gt, field3D, i, j)
    % determine the 2D cell index corresponding to the given (i,j) coordinates
    cell_ij = find(Gt.cells.ij(:,1) == i & Gt.cells.ij(:,2) == j);
    assert(numel(cell_ij) == 1, 'Expected exactly one cell for the given (i,j) coordinates');    

    % determining the corresponding 3D cells
    colcells3D = Gt.columns.cells(Gt.cells.columnPos(cell_ij):Gt.cells.columnPos(cell_ij+1)-1);

    % determining depth values of centroids
    z = Gt.parent.cells.centroids(colcells3D, 3);

    % extracting the field values at the corresponding 3D cells
    field_values = field3D(colcells3D);

    % plot the vertical profile
    figure;
    plot(field_values, z, '-o');
    set(gca, 'YDir', 'reverse'); % Reverse y-axis to have depth increasing
    xlabel('Field Value');
    ylabel('Depth (z)');
    title(sprintf('Vertical Profile of Field at (i=%d, j=%d)', i, j));
    grid on;
end
