%% ===================== Generic 3D Map (Any 3 Parameters + Isosurfaces + Slices + Switchable Colormap) =====================
% Dependency: GPRmodel_i, input order fixed as [Vr, MassRatio, Damp, Re]
load('GPR_Model1.mat'); 
%% ---------- User Settings ----------
% Axis variables (write in X/Y/Z order), and the remaining variable is fixed in "fixed"
axesNames = {'Vr','MassRatio','DampRatio'};   % e.g., {'Vr','MassRatio','Damp'} / {'Vr','Damp','Re'} etc.
fixed.Re  = 0.20;                        % Fixed value for the variable not used as an axis (example: fix Re)

% Value ranges for the three axes (denser = finer; number of prediction points = numel(vx)*numel(vy)*numel(vz))
Vals.Vr         = 3:0.15:14;
Vals.MassRatio  = 2:0.10:6;
Vals.DampRatio  = 0:0.002:0.04;
Vals.Re         = 0.10:0.02:0.30;        % Not used if this variable is fixed

% Isosurface control
isoVals      = 0.1:0.2:1.2;  % Cv isovalues to display
targetIso    = [];                      % Highlight target isovalue (empty [] means no highlight)
faceAlpha    = 0.75;                     % Transparency for normal isosurfaces
targetAlpha  = 0.70;                     % Transparency for target isosurface

% Optional orthogonal slices (to inspect internal volume data)
addSlices = false;                         % true/false
slicePos  = struct('X','mid','Y','mid','Z','mid');  % 'mid' or a numeric value (must be within axis ranges)

% Colormap modes:
% 'default'    → use MATLAB built-in colormap (specified by opts.cmap)
% 'smooth'     → use a "smooth gradient" colormap (blue-green, interpolated to 256 colors)
% 'customBase' → use your 9-color base and interpolate to 256 colors
opts.colorMode = 'customBase';
opts.cmap      = 'turbo';                % Only effective when colorMode='default'

% Colorbar range control (if not specified, use the isovalue range)
opts.caxisManual = false;
opts.caxisRange  = [0, 2];

% Figure window parameters
figPos = [80 80 500 300];

% Save output (optional)
saveFig.enable = false;                  % true/false
saveFig.path   = 'Cv_3D_isosurfaces.png';
saveFig.dpi    = 300;

%% ---------- Plot ----------
[Xv, Yv, Zv, Cv3d, info] = cv_volume_3d(axesNames, Vals, fixed, GPRmodel_i);
ax = plot_cv_isosurfaces(Xv, Yv, Zv, Cv3d, isoVals, targetIso, faceAlpha, targetAlpha, ...
                         addSlices, slicePos, figPos, info, opts);
% colormap(ax, cmap_base_bright(256, 1.30));

% camlight(ax,'headlight');
% camlight(ax,'right');
% camlight(ax,'left');
% lighting(ax,'gouraud');
material(ax,[0.3 0.85 0.5]);     % Balance specular + diffuse reflection

if saveFig.enable
    exportgraphics(ax.Parent, saveFig.path, 'Resolution', saveFig.dpi);
end
exportgraphics(gcf, 'Self_3D_MDVr.png', ...
'Resolution', 1200, ...
'BackgroundColor', 'none');

%% ===================== Function Section =====================
function [Xv, Yv, Zv, Cv3d, info] = cv_volume_3d(axesNames, Vals, fixed, model)
% Build a 3D grid and predict the Cv volume using the GPR model
% axesNames: 3-element cell array, e.g., {'Vr','MassRatio','Damp'}, order = X/Y/Z axes
% Vals     : vectors for each variable, field names match variable names
% fixed    : fixed value for the 4th variable not used as an axis, e.g., fixed.Re=0.2
% model    : GPRmodel_i, input order = [Vr, MassRatio, Damp, Re]

    varOrder = {'Vr','MassRatio','DampRatio','Re'};
    name2idx = containers.Map(varOrder,1:4);

    assert(numel(axesNames)==3, 'axesNames must contain 3 variable names.');
    assert(all(isKey(name2idx, axesNames)), 'Variable name must be one of Vr/MassRatio/Damp/Re.');

    % Find the variable name that is fixed
    fixedName = setdiff(varOrder, axesNames, 'stable');
    assert(numel(fixedName)==1, 'Exactly 1 variable must be fixed.');
    assert(isfield(fixed, fixedName{1}), 'Missing fixed value for the fixed variable in "fixed".');

    % Get axis vectors
    vx = Vals.(axesNames{1});
    vy = Vals.(axesNames{2});
    vz = Vals.(axesNames{3});

    % Generate 3D grid (ndgrid → (i,j,k) corresponds to (x,y,z))
    [Xv, Yv, Zv] = meshgrid(vx, vy, vz); 

    % Build query matrix (one 4D point per row: [Vr, MassRatio, Damp, Re])
    nP = numel(Xv);
    Xq = zeros(nP,4);

    % Fill the three axis variables into their corresponding columns
    Xq(:, name2idx(axesNames{1})) = Xv(:);
    Xq(:, name2idx(axesNames{2})) = Yv(:);
    Xq(:, name2idx(axesNames{3})) = Zv(:);

    % Fill the fixed variable column
    Xq(:, name2idx(fixedName{1})) = fixed.(fixedName{1});

    % Predict
    [Ypred,~,~] = predict(model, Xq);
    Cv3d = reshape(Ypred(:,1), size(Xv));

    % Return info
    info.axesNames = axesNames;         % Axis variable names (X,Y,Z)
    info.vx = vx; info.vy = vy; info.vz = vz;
    info.fixedName = fixedName{1};
    info.fixedVal  = fixed.(fixedName{1});
end


function ax = plot_cv_isosurfaces(Xv, Yv, Zv, Cv3d, isoVals, targetIso, faceAlpha, targetAlpha, ...
                                  addSlices, slicePos, figPos, info, opts)
% Plot multiple isosurfaces using isosurface; highlight target isosurface if provided.
% Optionally add orthogonal slices; colormap mode is switchable.

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on'); view(ax, 3);

    % ===== Select colormap mode =====
    switch opts.colorMode
        case 'default'
            colormap(ax, opts.cmap);

        case 'smooth'
            % Smooth gradient (blue-green example; you can replace with your preferred gradient)
            base2 = [ ...
                237,248,251;
                204,236,230;
                153,216,201;
                102,194,164;
                 65,174,118;
                 35,139, 69;
                  0, 90,  50] / 255;
            N2  = size(base2,1);
            xi2 = linspace(1, N2, 256);
            cmapSmooth = interp1(1:N2, base2, xi2, 'linear');
            colormap(ax, cmapSmooth);

        case 'customBase'
            % Your 9-color palette interpolated to 256 colors
            base = [ ...
                247,252,240;
                224,243,219;
                204,235,197;
                168,221,181;
                123,204,196;
                 78,179,211;
                 43,140,190;
                  8,104,172;
                  8, 64,129] / 255;
            N  = size(base,1);
            xi = linspace(1,N,256);
            cmapFromBase = interp1(1:N, base, xi, 'linear');
            colormap(ax, cmapFromBase);

        otherwise
            error('Unknown colormap mode opts.colorMode');
    end

    % Color axis limits
    if isempty(isoVals)
        clim(ax, [min(Cv3d(:)) max(Cv3d(:))]);
    else
        clim(ax, [min(isoVals) max(isoVals)]);
    end
    if isfield(opts,'caxisManual') && opts.caxisManual
        clim(ax, opts.caxisRange);
    end
    cb = colorbar(ax); 
    cb.Label.String = '$A^{*}_{CF}$';
    cb.FontSize = 10;           % Tick label font size
    cb.Label.FontSize = 12;     % Label font size
    cb.Label.Interpreter = 'latex';

    % Current colormap, used to sample colors based on iso values
    cmapCurrent = colormap(ax);
    if isempty(isoVals)
        isoVals = linspace(min(Cv3d(:)), max(Cv3d(:)), 5);
    end
    cmapIso = interp1( ...
        linspace(min(isoVals), max(isoVals), size(cmapCurrent,1)), ...
        cmapCurrent, ...
        isoVals, 'nearest', 'extrap');

    % ===== Draw isosurfaces =====
    for k = 1:numel(isoVals)
        val = isoVals(k);
        S = isosurface(Xv, Yv, Zv, Cv3d, val);
        if isempty(S.vertices), continue; end
        p = patch('Faces', S.faces, 'Vertices', S.vertices, ...
                  'FaceColor', cmapIso(k,:), ...
                  'EdgeColor', 'none', ...
                  'FaceAlpha', faceAlpha);
        isonormals(Xv, Yv, Zv, Cv3d, p); 
    end

    % ===== Highlight target isosurface (optional) =====
    if ~isempty(targetIso)
        S2 = isosurface(Xv, Yv, Zv, Cv3d, targetIso);
        if ~isempty(S2.vertices)
            p2 = patch('Faces', S2.faces, 'Vertices', S2.vertices, ...
                       'FaceColor', [0.85 0.30 0.30], ...  % Deep red highlight
                       'EdgeColor', 'none', ...
                       'FaceAlpha', targetAlpha);
            isonormals(Xv, Yv, Zv, Cv3d, p2);
        end
    end

    % ===== Orthogonal slices (optional) =====
    if addSlices
        xs = pickSlicePos(info.vx, slicePos.X);
        ys = pickSlicePos(info.vy, slicePos.Y);
        zs = pickSlicePos(info.vz, slicePos.Z);
        hs = slice(ax, Xv, Yv, Zv, Cv3d, xs, ys, zs);
        set(hs, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end

    % View/lighting/axes
    axis(ax,'tight'); 
    axis(ax,'normal');        % ✅ Key: show each dimension in its own scale!
    
    set(ax, 'DataAspectRatioMode','auto');
    set(ax, 'PlotBoxAspectRatioMode','auto');
    
    view(ax, 3);
    grid(ax,'on'); box(ax,'on');
    camlight(ax,'headlight'); lighting(ax,'gouraud');
    ax.FontSize = 10;
    xlabel(ax, info.axesNames{1}, 'FontSize', 12);
    ylabel(ax, info.axesNames{2}, 'FontSize', 12);
    zlabel(ax, info.axesNames{3}, 'FontSize', 12);

    % title(ax, sprintf('%s=%.4d', info.fixedName, info.fixedVal*100000), ...
    %       'FontSize', 14);
end


function v = pickSlicePos(vec, spec)
% spec='mid' selects the midpoint; or a numeric value (clamped to the vector range)
    if (ischar(spec) || isstring(spec)) && strcmpi(spec,'mid')
        v = median(vec);
    else
        v = min(max(double(spec), min(vec)), max(vec));
    end
end

function cmapBright = cmap_base_bright(levels, brightenFactor)
    if nargin<1, levels = 256; end
    if nargin<2, brightenFactor = 1.25; end     %  Brightness boost factor

    base = [ ...
        247,252,240;
        224,243,219;
        204,235,197;
        168,221,181;
        123,204,196;
         78,179,211;
         43,140,190;
          8,104,172;
          8, 64,129] / 255;

    % Brighten: multiply the color matrix by a factor > 1
    baseBright = min(base * brightenFactor, 1);  %  Avoid exceeding 1

    xi = linspace(1, size(baseBright,1), levels);
    cmapBright = interp1(1:size(baseBright,1), baseBright, xi, 'linear');
end
