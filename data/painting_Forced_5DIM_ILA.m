% Define four X2 values (IL-direction amplitude A*il)
X2_values = 0.1:0.1:0.3;
Re = 0.1:0.01:0.3;
f_fine = (0.05:0.001:0.35)';
X1 = 0.1:0.1:1.5;    % amplitude range
Y1 = 0.01:0.01:0.35; % frequency range
THETA = 320;         % fixed phase angle at 320°

% Create figure window
figure;
set(gcf, 'Position', [100, 100, 900, 220], 'Color', 'w');
set(groot,'DefaultAxesFontName','Segoe UI');
set(groot,'DefaultTextFontName','Segoe UI');
% Create a black-to-red gradient color scheme
colors = [0.00, 0.00, 0.00   % pure black
          1.0, 0.7, 0.7      % light red
          0.70, 0.00, 0.00   % dark red
          0.70, 0.00, 0.00   % deep red
          0.77, 0.12, 0.23];   % red
positions = [0, 0.03, 0.6, 0.85, 1.0];

cmap = interp1(positions, colors, linspace(0,1,256));
colormap(cmap); % set global colormap

% Precompute max A value for unified axis limits
max_A = 0;
min_A = inf;
all_A_data = cell(1, length(X2_values));

for x2_idx = 1:length(X2_values)
    current_X2 = X2_values(x2_idx);
    A_at_ce0_plot = zeros(length(f_fine), length(Re));
    
    for m = 1:length(Re)
        [x1,y1,x2,theta,re] = ndgrid(X1,Y1,current_X2,THETA,Re(m));
        Xtrain = [x1(:),y1(:),x2(:),theta(:),re(:)];

        load('GPR_Model1.mat')
        [Ytrain, ~, ~] = predict(GPRmodel_i, Xtrain);
        K = Ytrain(:,1);

        % Reconstruct coefficient matrix
        Cv = reshape(K, [length(X1), length(Y1)]);
        ce = Cv;
        
        %% Zero-crossing detection
        A_at_ce0_raw = zeros(1, length(Y1));
        
        for j = 1:length(Y1)
            ce_col = ce(:,j);
            cross_points = [];
            
            for i = 1:(length(X1)-1)
                v1 = ce_col(i);
                v2 = ce_col(i+1);
                
                if v1*v2 <= 0 && ~(v1==0 && v2==0)
                    A_interp = X1(i) + (0 - v1)*(X1(i+1)-X1(i))/(v2 - v1);
                    
                    if A_interp >= X1(1) && A_interp <= X1(end)
                        cross_points(end+1) = A_interp;
                    elseif A_interp > X1(end) && v2 < 0
                        cross_points(end+1) = A_interp;
                    end
                end
            end
            
            if ~isempty(cross_points)
                A_at_ce0_raw(j) = max(cross_points);
            else        
                A_at_ce0_raw(j) = 0;
            end
            if Y1(j)<0.06
                A_at_ce0_raw(j) = 0;
            end
        end
        
        A_at_ce0_fine = interp1(Y1, A_at_ce0_raw, f_fine, 'pchip', 'extrap');
        A_at_ce0_fine(A_at_ce0_fine < 0) = 0;
        A_at_ce0_plot(:,m) = A_at_ce0_fine;
    end
    
    % Save data for unified axis limits
    all_A_data{x2_idx} = A_at_ce0_plot;
    current_max = max(A_at_ce0_plot(:));
    current_min = min(A_at_ce0_plot(:));
    
    if current_max > max_A
        max_A = current_max;
    end
    if current_min < min_A
        min_A = current_min;
    end
end

% Compute unified axis limits
z_lim = [min_A * 0.95, max_A * 1.05];
if z_lim(1) < 0
    z_lim(1) = 0; % ensure amplitude is not less than 0
end

% Create all subplots
for x2_idx = 1:length(X2_values)
    current_X2 = X2_values(x2_idx);
    A_at_ce0_plot = all_A_data{x2_idx};
    
    %% Create subplot
    subplot(1, 3, x2_idx);
    
    % Generate grid data
    [F_mesh, R_mesh] = meshgrid(f_fine, Re*1e5);
    A_mesh = A_at_ce0_plot';
    
    % Create surface plot - use global colormap
    surf(F_mesh, R_mesh, A_mesh, 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud','FaceAlpha', 0.6);
    shading interp
    
    % Set unified color axis range
    clim(z_lim);
    
    % Set unified view and axis limits
    view(-40, 28)
    xlim([min(f_fine), max(f_fine)])
    ylim([min(Re)*1e5, max(Re)*1e5])
    zlim(z_lim) % unified z-axis range
    
    % Add title and axis labels
    title(['$A^{*}_{IL} = $', num2str(current_X2)], 'FontSize', 12, 'Interpreter','latex')
    % title(['\A*_IL = ', num2str(current_X2)], 'FontSize', 12)
    if x2_idx == 1
        zlabel('$A^{*}_{CF}$','FontSize', 12, 'Interpreter','latex')
        ylabel('Re', 'FontSize', 12)
    else
        set(gca, 'YTickLabel', [])
    end
    
    % Simplify axis annotations
    set(gca, 'FontSize', 10)
    xlabel('$f^{*}$', 'FontSize', 12, 'Interpreter','latex')
    
    grid on
    box on
end

exportgraphics(gcf, 'Forced_5D_AIL.png', ...
'Resolution', 1800, ...
'BackgroundColor', 'none');

%% Add shared colorbar (move down)
% % Compute colorbar position (bottom of the entire figure)
% cbar_height = 0.03; % colorbar height
% cbar_bottom = 0.05; % distance from bottom
% cbar_width = 0.7;   % colorbar width
% cbar_left = (1 - cbar_width)/2; % horizontally centered
% 
% % Create global colorbar
% c = colorbar('Position', [cbar_left, cbar_bottom, cbar_width, cbar_height], ...
%              'Orientation', 'horizontal');
% c.Label.String = 'Critical amplitude A* (C_v=0)';
% c.Label.FontSize = 12;
% c.Label.FontWeight = 'bold';
% clim(z_lim) % set global color axis range
% 
% % Adjust subplot positions to leave space for the colorbar
% subplot_margin = 0.15; % leave space at bottom of subplots
% for i = 1:4
%     h = subplot(1,4,i);
%     pos = get(h, 'Position');
%     set(h, 'Position', [pos(1), pos(2)+subplot_margin/2, pos(3), pos(4)-subplot_margin]);
% end

% Add overall title
% sgtitle(['Critical amplitude distribution at phase angle \theta = ', num2str(THETA), '°'], ...
%         'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.3]);
