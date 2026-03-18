function h = plot_network(p,e,t,e_coarse,t_coarse)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

gcf; hold on
h = [];

% Plot fine interface mesh (fine edges)
for i = 1:length(e(1,:))
    x = [p(1,e(1,i)), p(1,e(2,i))];
    y = [p(2,e(1,i)), p(2,e(2,i))];
    h_temp = plot(x,y,'-b','LineWidth',2);%,'MarkerSize',3);
    h = [h; h_temp];
end

if nargin >= 3
    % Plot bulk (fine) mesh
    for i = 1:length(t(1,:))
        x = [p(1,t(1,i)) p(1,t(2,i)) p(1,t(3,i)) p(1,t(1,i))];
        y = [p(2,t(1,i)) p(2,t(2,i)) p(2,t(3,i)) p(2,t(1,i))];
        h_temp = plot(x,y,'r');
        h = [h; h_temp];
    end
end


% Plot interfaces
for i = 1:length(e(1,:))
    x = [p(1,e(1,i)) p(1,e(2,i))];
    y = [p(2,e(1,i)) p(2,e(2,i))];
    h_temp = plot(x,y,'b');
    h = [h; h_temp];
end

    % % Plot coarse interface mesh
    % for i = 1:length(e_coarse(1,:))
    %     x = [p(1,e_coarse(1,i)), p(1,e_coarse(2,i))];
    %     y = [p(2,e_coarse(1,i)), p(2,e_coarse(2,i))];
    %     h_temp = plot(x,y,'-g','LineWidth',10);
    %     h = [h; h_temp];
    % end

if nargin == 5 % Coarse mesh should also be plotted
    % % Plot coarse bulk mesh
    % for i = 1:length(t_coarse(1,:))
    %     x = [p(1,t_coarse(1,i)) p(1,t_coarse(2,i)) p(1,t_coarse(3,i))];% p(1,t_coarse(1,i))];%
    %     y = [p(2,t_coarse(1,i)) p(2,t_coarse(2,i)) p(2,t_coarse(3,i))];% p(2,t_coarse(1,i))];
    %     h_temp = plot(x,y,'r','LineWidth',2);
    %     h = [h; h_temp];
    % end

    % Plot coarse bulk mesh
    % Do not plot the ones interfering with interfacess
    for i = 1:length(t_coarse(1,:))
        if ~any(all([p(4,t_coarse(1,i)); p(4,t_coarse(2,i))] == e_coarse)) && ~any(all([p(4,t_coarse(2,i)); p(4,t_coarse(1,i))] == e_coarse))
            x = [p(1,t_coarse(1,i)) p(1,t_coarse(2,i))];% p(1,t_coarse(1,i))];%
            y = [p(2,t_coarse(1,i)) p(2,t_coarse(2,i))];% p(2,t_coarse(1,i))];
            h_temp = plot(x,y,'r','LineWidth',2);
            h = [h; h_temp];
        end
        if ~any(all([p(4,t_coarse(3,i)); p(4,t_coarse(2,i))] == e_coarse)) && ~any(all([p(4,t_coarse(2,i)); p(4,t_coarse(3,i))] == e_coarse))
            if i == 36, disp('36'),end
            x = [p(1,t_coarse(2,i)) p(1,t_coarse(3,i))];% p(1,t_coarse(1,i))];%
            y = [p(2,t_coarse(2,i)) p(2,t_coarse(3,i))];% p(2,t_coarse(1,i))];
            h_temp = plot(x,y,'r','LineWidth',2);
            h = [h; h_temp];
        end
        if ~any(all([p(4,t_coarse(1,i)); p(4,t_coarse(3,i))] == e_coarse)) && ~any(all([p(4,t_coarse(3,i)); p(4,t_coarse(1,i))] == e_coarse))
            x = [p(1,t_coarse(1,i)) p(1,t_coarse(3,i))];% p(1,t_coarse(1,i))];%
            y = [p(2,t_coarse(1,i)) p(2,t_coarse(3,i))];% p(2,t_coarse(1,i))];
            h_temp = plot(x,y,'r','LineWidth',2);
            h = [h; h_temp];
        end
    end
    % % Plot coarse interface mesh
    % for i = 1:length(e_coarse(1,:))
    %     x = [p(1,e_coarse(1,i)), p(1,e_coarse(2,i))];
    %     y = [p(2,e_coarse(1,i)), p(2,e_coarse(2,i))];
    %     h_temp = plot(x,y,'b','MarkerSize',5);
    %     h = [h; h_temp];
    % end
end

% axis([0 1 0 1])
% xticks([0 1])
% yticks([0 1])
axis equal
axis off

end






