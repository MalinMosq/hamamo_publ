% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

%% Generate domain
clc,clear
if ~exist('utils/data_domain', 'dir'), mkdir('utils/data_domain'), end
if ~exist('utils/data_mesh', 'dir'), mkdir('utils/data_mesh'), end
addpath utils/

% Test case 4 lines: nlines = 7, rng_no = 33
% Test case 1 line: nlines = 1, rng_no = 10
% Test case 8 lines: nlines = 10, rng_no = 26
n_lines = 7;
finite = 0; % corresponds to graphtype = plp (0, infinite lines = convex subdomains) or fl (1, finite lines = nonconvex subdomains)
tol = 1e-8;
plot_domain = 1; % Slow when n is high!
rng_no = 33;

disp('------------------ Generate domain ------------------')
generate_domain(n_lines,finite,tol,plot_domain,rng_no)
% Saves to file:
% - edges: which points are connected by edges
% - points: x- and y-coordinates of the points in the domain
% - subdomains_e: describes all the subdomains referring to the edges
% - subdomains_p: describes all the subdomains referring to the points

 %% Generate mesh (python)
% disp('------------------- Generate mesh -------------------')
% % Check compatibility with your version of MATLAB and Python here:
% % https://se.mathworks.com/support/requirements/python-compatibility.html
% % Link to the executable on your computer:
% disp(pyenv)
% % Python 3.13 is not yet supported by Matlab. Run py-file from terminal
% % instead. Python 3.9 is installed in Matlab but gmsh is not installed
% % for this version of Python.
% pyenv('Version','/Library/Frameworks/Python.framework/Versions/3.13/bin/python3');
% % pyenv('ExecutionMode','OutOfProcess');
% 
% pyrunfile('utils/generate_mesh.py')
% % Alternative: pause here and run python file manually
% % Saves to file:
% % - mesh1.m, containing a coarse mesh

%% Refine mesh
disp('-------------------- Refine mesh --------------------')
h = 2;
plot_meshrefinement = 0;
for H = 0:h-2
    H
    refine_mesh(h,H, plot_meshrefinement)
end
% Saves to file:
% - mesh2.m, mesh4.m, etc depending on nbr_meshes, each file a refinement
% of the previous one


%% Pose problem on geometry and solve it
clc
disp('--------------- Pose problem and solve ---------------')
addpath utils/
tol = 1e-8;
h = 5; H = h-2;
pde_problem(h,H,tol);








%% Below are some scripts for producing figures etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decay of basis functions
clear

load('lod_basis_h6_H2_B1_ell4.mat')
whos
Phi_enh_loc = Phi;
n_tri = size(t_coarse,2); n_edg = size(e_coarse,2);

figure('color','white')
set(gcf,'Position',[0,0,1500,1500])
% subplot(1,2,1)
% h = plot_network(p,e_coarse,t_coarse);
%zdatam(h,-30)

for i = 1:size(Phi_enh_loc,2)
    if sum(Phi_enh_loc(:,i)) < 1e-14
        continue
    end
    
h = plot_network(p,e_coarse,t_coarse);

    x = reshape(p(1,t(1:3,:)),3,[]);
    y = reshape(p(2,t(1:3,:)),3,[]);
    z = abs(reshape(Phi_enh_loc(t(1:3,:),i),3,[]));
    % h1 = patch(x,y,z,z,'facecolor','interp','edgecolor','none');%,'facealpha',.5);

    % Plot bulk parts
    h1 = trisurf(t(1:3,:)',p(1,:),p(2,:),log(abs(Phi_enh_loc(:,i))),'facecolor','interp','edgecolor','none');
    % Plot interfaces
    h2 = [];
    for j = 1:size(e,2)
        h_temp = patch([p(1,e(1,j)) p(1,e(2,j))],...
            [p(2,e(1,j)) p(2,e(2,j))],...
            ones(size([p(2,e(1,j)) p(2,e(2,j))])),...%log([abs(Phi_enh_loc(e(1,j),i)) abs(Phi_enh_loc(e(2,j),i))])
            log([abs(Phi_enh_loc(e(1,j),i)) abs(Phi_enh_loc(e(2,j),i))]),...
            'FaceColor','none','EdgeColor','interp','linewidth',2);
        h2 = [h2, h_temp];
    end % for
    axis([0 1 0 1])
    cbh = colorbar; cbh.FontSize = 16;
    ticks = get(cbh,'Ticks');
    ticks = unique(floor(log10(exp(ticks))));
    ticks = setdiff(log(10.^ticks),0);
    ticks = [ticks, 0];
    set(cbh,'Ticks',ticks)
    ticks = exp(ticks);
    set(cbh,'TickLabels',(ticks))
disp('------')

% subplot(1,2,2)
% x = p(1,:); y = p(2,:);
% z = (abs(Phi_enh_loc(:,i)))/max(abs(Phi_enh_loc(:,i))); % c = (abs(Phi_enh_loc(:,i)))/max(abs(Phi_enh_loc(:,i)))
% z(z==0) = NaN;
% p_inter_logical = zeros(size(p,2),1); p_inter_logical(unique(e(1:2,:))) = 1;
% bulk = ~p_inter_logical;
% % z = sin(x+y);
% % max(z)
% trisurf(t(1:3,:)',x,y,log(z),'edgecolor','none') %,'facecolor','interp' %(bulk)
%     h2 = [];
%     for j = 1:size(e,2)
%         h_temp = patch([p(1,e(1,j)) p(1,e(2,j))],...
%             [p(2,e(1,j)) p(2,e(2,j))],...
%             log([abs(Phi_enh_loc(e(1,j),i)) abs(Phi_enh_loc(e(2,j),i))]),...
%             log([abs(Phi_enh_loc(e(1,j),i)) abs(Phi_enh_loc(e(2,j),i))]),...
%             'FaceColor','none','EdgeColor','interp','linewidth',2);
%         h2 = [h2, h_temp];
%     end % for


    a = axes;
    %t1 = title(['Basis ', num2str(i)]);
    a.Visible = 'off'; % set(a,'Visible','off');
    %t1.Visible = 'on'; % set(t1,'Visible','on');
    view(2)


    % Image = getframe(gcf);
    % imwrite(Image.cdata, ['figures/bases_h6_H2/LOD_basis_',num2str(i),'.jpg'],'quality',100);
    % imwrite(Image.cdata, ['figures/bases_h6_H2/LOD_basis_',num2str(i),'.png']);
    exportgraphics(gcf,['figures/lod_basis_h6_H2_B1_ell4/LOD_basis_',num2str(i),'.png'])
    % matlab2tikz(['figures/bases_h6_H2/LOD_basis_',num2str(i),'.tex'])bases_h6_H2_B100_ell4
    % pause(1)
    % delete(h1),delete(h2),delete(h4),delete(h5)
    delete(h1), delete(h2)
    delete(h)
end



%% 

load('sol_new_complicated1.1_h5.mat')

size(A)
size(u_fem)

h = plot_solution(p,e_inter,t,u_fem);

exportgraphics(gca,'figures/sol_new_complicated1.1_h6.png')

%% Inzoomad bild
clc
figure(),hold on
h = 6;
for H = 2%0:4
    [p,e,t,p_inter_logical,p_boundary_logical,e_inter,e_coarse,t_coarse] = ...
        read_variables(h,H,tol);
    

    handle = [];
    xmin = .8; xmax = .95; ymin = .45; ymax = .6; xmid = .865; ymid = .515;
    
    for i = 1:length(t(1,:))
        x = [p(1,t(1,i)) p(1,t(2,i)) p(1,t(3,i)) p(1,t(1,i))];
        y = [p(2,t(1,i)) p(2,t(2,i)) p(2,t(3,i)) p(2,t(1,i))];
        % h_temp = plot(x,y,'r');
        % handle = [handle; h_temp];

        %if and(all(and(x >= xmin, x <= xmax)),all(and(y >= ymin, y <= ymax)))
        if all( ((x-xmid).^2 + (y - ymid).^2) <= .06^2 )
            plot(x-xmid,y-ymid,'r')
        end
    end
    % subplot(3,2,H+1), hold on
    % title(num2str(H))
    % Plot coarse bulk mesh
    for i = 1:length(t_coarse(1,:))
        x = [p(1,t_coarse(1,i)) p(1,t_coarse(2,i)) p(1,t_coarse(3,i)) p(1,t_coarse(1,i))];
        y = [p(2,t_coarse(1,i)) p(2,t_coarse(2,i)) p(2,t_coarse(3,i)) p(2,t_coarse(1,i))];
        % h_temp = plot(x,y,'r','LineWidth',1);
        % handle = [handle; h_temp];
        if and(any(and(x >= xmin, x <= xmax)),any(and(y >= ymin, y <= ymax)))
            plot(x-xmid,y-ymid,'r','LineWidth',2)
        end
    end
    % % Plot coarse interface mesh
    % for i = 1:length(e_coarse(1,:))
    %     x = [p(1,e_coarse(1,i)), p(1,e_coarse(2,i))];
    %     y = [p(2,e_coarse(1,i)), p(2,e_coarse(2,i))];
    %     h_temp = plot(x,y,'ob','MarkerSize',5);
    %     handle = [handle; h_temp];
    % end
    % % Plot interfaces
    % for i = 1:length(e(1,:))
    %     x = [p(1,e(1,i)) p(1,e(2,i))];
    %     y = [p(2,e(1,i)) p(2,e(2,i))];
    %     h_temp = plot(x,y,'b', 'LineWidth',2);
    %     handle = [handle; h_temp];
    % end
end
axis equal, axis off
% text(0,-.05,'0','FontSize',12)
% text(1,-.05,'1','FontSize',12)
% text(-.05,1,'1','FontSize',12)
% text(-.05,0,'0','FontSize',12)

%matlab2tikz('Domain_4lines.tex')

%% Plot mesh
h = 5;
H = 1;
[p,e,t,p_inter_logical,p_boundary_logical,e_inter,e_coarse,t_coarse] = ...
    read_variables(h,H,tol);
plot_network(p,e_coarse,t_coarse);
exportgraphics(gca,'figures/domain_complicated1.1.png')

