function [G, rock, Gt, rockVE, transMult, G_orig, rock_orig] = setup_testgrid(varargin)

    opt = merge_options(struct('uniform', false), varargin{:});
    
    % load grid slice and rock properties:
    ld = load('grids');
    Gs = ld.Gs;                % The full 3D grid slice
    Gmid = ld.Gmid;            % The middle section of the slice, sandwiched
                               % between two low-permeability layers
    rock = ld.rocks;           % Rock properties for the full grid slice
    rock = ld.rocks_mid;   % Rock properties for the middle section of the slice

    if opt.uniform
        rock.perm(:) = 100; % uniform permeability of 100 mD
        rock.permz(:) = 100; % uniform vertical permeability of 10 mD
        rock.poro(:) = 0.15; % uniform porosity of 20%
        rock.ntg(:) = 1; % uniform net-to-gross of 1 (fully net)
        rock.facies(:) = 1; % uniform facies of 1 (arbitrary choice)
    end
    
    % convert units
    rock.perm = rock.perm * milli * darcy; % from mD to Darcy
    rock.permz = rock.permz * milli * darcy; % from mD to Darcy
    rock.poro = rock.poro; % porosity is unitless

    % removing rightmost column of Gmid, which has practically zero permeability,
    % thus rendering the use of open boundary conditions meaningless.
    tmp = zeros(Gmid.cartDims);
    tmp(Gmid.cells.indexMap) = 1:Gmid.cells.num;
    rightmost_col_cells = unique(tmp(end,:,:));
    rightmost_col_cells = rightmost_col_cells(rightmost_col_cells > 0);
    keep_cells = setdiff(1:Gmid.cells.num, rightmost_col_cells);
    
    Gmid = extractSubgrid(Gmid, keep_cells);

    rock.perm = rock.perm(keep_cells);
    rock.permz = rock.permz(keep_cells);
    rock.poro = rock.poro(keep_cells);
    rock.ntg = rock.ntg(keep_cells);
    rock.facies = rock.facies(keep_cells);
    
    % Define a simple vertical equilibrium grid
    [Gt, G, transMult, discarded_cells] = topSurfaceGrid(Gmid, 'discard_below_holes', false);
    
    G_orig = Gmid; % store original grid before discarding cells
    rock_orig = rock; % store original rock properties before discarding cells
    
    % loop over fields in rocks_mid and remove values corresponding to discarded
    % cells
    rock_fields = fieldnames(rock);
    for k = 1:numel(rock_fields)
        field = rock_fields{k};
        rock.(field)(discarded_cells) = [];
    end
    
    rockVE = averageRock(rock, Gt);
end
