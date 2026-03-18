function generate_domain(n,finite,tol,plotta,rng_no)
% M. Hauck, A. Målqvist, M. Mosquera, 2026
% A LOCALIZED ORTHOGONAL DECOMPOSITION METHOD FOR HETEROGENEOUS MIXED-DIMENSIONAL PROBLEMS

if nargin == 5
    rng(rng_no)
end

if finite
    len = .2;
    [p,e,conv] = finite_lines(n,tol,len);
else
    [p,e,conv] = infinite_lines(n,tol);
end

[p,e] = remove_duplicates(p,e);
remove_duplicates_time = toc;


if plotta
    plot_network(p,e); % plot of network
end

tic
[cycles_p, hanging_edges] = findshortestcyclebasis(p,e,tol,conv);
findcycles_time = toc;

if ~isempty(hanging_edges)
    disp('OBS hanging edges that are not part of a subdomain:')
    disp(hanging_edges)
    if plotta
        hanging_edges = [hanging_edges(:,1);hanging_edges(:,2)];
    end
end

cycles_e = cycles_points2edges(e,cycles_p);
cycles_points2edges_time = toc;

save_variables(p,e,cycles_p, cycles_e)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,e,conv] = finite_lines(n,tol,len)

[xstart,xstop,ystart,ystop,mids] = fl(n,len);
[p,e] = discretizefl(xstart,xstop,ystart,ystop,mids,len);
p(:,1)

p(:,end+1:end+4) = [0 0 1 1; 0 1 0 1];

% Add x = 1
p_ind = find(abs(p(1,:) - 1) < tol);
[~,I] = sort(p(2,p_ind));
p_ind = p_ind(I);
e = [e, [p_ind(1:end-1);p_ind(2:end)]];
% Add x = 0
p_ind = find(abs(p(1,:)) < tol);
[~,I] = sort(p(2,p_ind));
p_ind = p_ind(I);
e = [e, [p_ind(1:end-1);p_ind(2:end)]];
% Add y = 1
p_ind = find(abs(p(2,:) - 1) < tol);
[~,I] = sort(p(1,p_ind));
p_ind = p_ind(I);
e = [e, [p_ind(1:end-1);p_ind(2:end)]];
% Add y = 0
p_ind = find(abs(p(2,:)) < tol);
[~,I] = sort(p(1,p_ind));
p_ind = p_ind(I);
e = [e, [p_ind(1:end-1);p_ind(2:end)]];

conv = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,e,conv] = infinite_lines(n,tol)
%n is number of randomly distributed infinity lines in the unit square
%N is the number of coarse elements of the mesh T_H is one coordinate
global display_times

tic
[x1,x2,y1,y2]=plp(n,tol);%n random infinite lines in the unit square
x1(abs(x1 - 1) < tol) = 1; x1(abs(x1) < tol) = 0;
x2(abs(x2 - 1) < tol) = 1; x2(abs(x2) < tol) = 0;
y1(abs(y1 - 1) < tol) = 1; y1(abs(y1) < tol) = 0;
y2(abs(y2 - 1) < tol) = 1; y2(abs(y2) < tol) = 0;
%%%%%%%% Nedan är för test case %%%%%%%%%
if and(n == 7, rng().Seed == 33)
    x1(1) = []; x2(1) = []; y1(1) = []; y2(1) = [];
elseif and(n == 1, rng().Seed == 1)
    x1 = 0; x2 = 1; y1 = .5; y2 = .5;
elseif and(n == 1, rng().Seed == 3)
    x1 = .3; x2 = .3; y1 = 0; y2 = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [x1; 0; 0; 0; 1]; x2 = [x2; 1; 1; 0; 1]; % Add boundary
y1 = [y1; 0; 1; 0; 0]; y2 = [y2; 0; 1; 1; 1]; % Add boundary
[p,e]=discretize(x1,x2,y1,y2,tol); %intersections are nodes (p) connected by edges (e)
generation_time = toc;
if display_times
    disp(['generation_time = ', num2str(generation_time)])
end
conv = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cycles_e = cycles_points2edges(e,cycles_p)
max_length = size(cycles_p,2)-1;
cycles_n = size(cycles_p,1);
cycles_e = zeros(cycles_n,max_length);

for i = 1:cycles_n
    cycle_length = sum(cycles_p(i,:)>0) - 1;
    edge_cycle = zeros(1,cycle_length);
    for j = 1:cycle_length
        edge_cycle(j) = find(all(e == sort(cycles_p(i,j:j+1))',1));
    end
    cycles_e(i,1:cycle_length) = edge_cycle;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_variables(p,e,cycles_p, cycles_e)

% Save all the points
fileID = fopen('utils/data_domain/points.txt','w');
for i = 1:length(p(1,:))
    fprintf(fileID, '%.8f %.8f\n',p(:,i)');
end
fclose(fileID);

% Save all the edges
e = uint16(e);
fileID = fopen('utils/data_domain/edges.txt','w');
for i = 1:length(e(1,:))
    fprintf(fileID, '%u %u\n',e(:,i)');
end
fclose(fileID);

% Save all the subdomains (cycles)
fileID = fopen('utils/data_domain/subdomains_e.txt','w');
for i = 1:size(cycles_e,1)
    for j = 1:size(cycles_e,2)
        fprintf(fileID, '%u ',cycles_e(i,j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

% Save all the subdomains (cycles)
fileID = fopen('utils/data_domain/subdomains_p.txt','w');
for i = 1:size(cycles_p,1)
    for j = 1:size(cycles_p,2)
        fprintf(fileID, '%u ',cycles_p(i,j));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,e] = remove_duplicates(p,e)

% Find duplicates and their first occurrence
p_round = round(p,8);
val = [];
for i = 1:size(p,2)
    occurrences = all(p_round == p_round(:,i),1);
    if sum(occurrences) > 1 && sum(occurrences(1:i)) ~= 1
        val = [val; [find(occurrences,1), i]];
    end
end

% Change point indices in e so that only the correct points are there
if ~isempty(val)
    e = changem(e,val(:,1),val(:,2));
end

% Remove edges going to the same point
e(:,logical(e(1,:) == e(2,:))) = [];

% Change point indices in e so that they refer to right point in p
% and remove superfluous points in p
if ~isempty(val)
    removed_points = unique(val(:,2));
    for i = 1:length(removed_points)
        e(e > removed_points(end+1-i)) = e(e > removed_points(end+1-i)) - 1;
    end
    p(:,removed_points) = [];
end

% Remove duplicates in e (if there are any after merging points)
e = sort(e,1);
e = unique(e','rows')';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cycles, hanging_edges] = findshortestcyclebasis(p,e,tol,conv)

cycles = [];
cycles_ordered = [];
untreated_edges = e;
hanging_edges = [];

for kk = 1:size(e,2)
    % e(:,1) is the edge going from e(1,1) to e(2,1)
    % disp('---------')
    p1 = untreated_edges(1,1);
    p2 = untreated_edges(2,1);
    % Find connected points
    [rows,cols] = find(e == p2);

    rows = mod(rows,2)+1; % Choose connected points
    p_next = e(sub2ind(size(e),rows,cols));
    p_next = p_next(p_next ~= p1);
    % Find angles
    alpha_start = findangles(p,[p1 p2],p_next);
    [~, minindex] = min(alpha_start);
    loop_min = [p1 p2 p_next(minindex)];
    [~, maxindex] = max(alpha_start);
    loop_max = [p1 p2 p_next(maxindex)];
    if isempty(p_next) % In large networks, some edges can disappear
        untreated_edges(:,1) = [];
        hanging_edges = [hanging_edges; p1, p2];
        if isempty(untreated_edges), break, end
        continue
    end
    % Break if we're building on an already existing cycle
    % Convex areas => we'll notice this directly after finding p3
    if conv
        if any(all([any(cycles==loop_min(1),2),...
                any(cycles==loop_min(2),2),...
                any(cycles==loop_min(3),2)],2))
            alpha_start(minindex) = pi;
        end
        if any(all([any(cycles==loop_max(1),2),...
                any(cycles==loop_max(2),2),...
                any(cycles==loop_max(3),2)],2))
            alpha_start(maxindex) = pi;
        end
    end
    if alpha_start(minindex) < (pi - tol)
        % Always choose the minimum angle!
        [cycles, cycles_ordered] = loop_builder(cycles,cycles_ordered,loop_min,p,e,conv,'min');
    end
    if alpha_start(maxindex) > (pi + tol)
        % Always choose the maximum angle!
        [cycles, cycles_ordered] = loop_builder(cycles,cycles_ordered,loop_max,p,e,conv,'max');
    end
    untreated_edges(:,1) = [];
    if isempty(untreated_edges), break, end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cycles, cycles_ordered] = loop_builder(cycles,cycles_ordered,loop,p,e,conv,maxmin)

for j = 1:length(p)
    % Find connected points
    [rows,cols] = find(e == loop(end));
    rows = mod(rows,2)+1; % Choose connected points
    p_next = e(sub2ind(size(e),rows,cols));
    p_next = p_next(p_next ~= loop(end-1));

    % Find angles
    alpha = findangles(p,loop,p_next);
    eval(['[~, maxindex] = ' maxmin '(alpha);'])
    loop = [loop p_next(maxindex)];
    % Break if the loop is finished
    if loop(end) == loop(1)
        if size(cycles,2) < length(loop)
            cycles = [cycles, ...
                zeros(size(cycles,1), length(loop)-size(cycles,2))];
            if ~conv
                cycles_ordered = [cycles_ordered, ...
                    zeros(size(cycles_ordered,1), ...
                    length(loop)-size(cycles_ordered,2))];
            end
        end
        loop = [loop, ...
            zeros(1,size(cycles,2) - length(loop))];
        if ~conv
            loop_ordered = [flip(unique(loop)),...
                zeros(1,size(cycles_ordered,2) - length(unique(loop)))];
            if any(all(cycles_ordered == loop_ordered, 2))
                % Cycle already exists in cycle variable
                break
            end
            cycles_ordered(end+1,:) = loop_ordered;
        end
        cycles(end+1,:) = loop;
        break
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha = findangles(p,loop,p_next)
% As per
% https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
alpha = pi*ones(size(p_next));
for i = 1:length(p_next)
    p3 = p_next(i);
    v_old = p(:,loop(end-1)) - p(:,loop(end));
    v_new = p(:,p3) - p(:,loop(end));
    dot_v = v_old'*v_new;
    det_v = det([v_old v_new]);
    alpha(i) = atan2(det_v,dot_v);
end
alpha = mod(alpha, 2*pi);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,e] = discretize(x1,x2,y1,y2,tol)
% x1|y1 and x2|y2 are start and stop coord for lines in a unit square 
% all lines pass through interior points (on frame allowed)
% p is coords for all intersections and e are the edges of the network

%p matrix
n=length(x1);N=nan(n,n+2);p=[];
for i=1:n
    p=[p [x1(i);y1(i)]];
    N(i,1)=size(p,2);
    %old nodes
    for j=1:i-1
        N(i,j+1)=N(j,i+1);
    end
%     if i == 57
%         disp([x1(i) y1(i) x2(i) y2(i)])
%         disp('---')
%     end

    for j=i+1:n
        %compute angles
        thetai=atan((y2(i)-y1(i))./(x2(i)-x1(i)));
        thetaj=atan((y2(j)-y1(j))./(x2(j)-x1(j)));
        %coord where/if line i and j cut
        R=[cos(thetaj) -cos(thetai);sin(thetaj) -sin(thetai)];


        if abs(det(R))>1e-11 %not parallell
            s=R\[x1(i)-x1(j);y1(i)-y1(j)];
            pxij=x1(i)+s(2)*cos(thetai);
            pyij=y1(i)+s(2)*sin(thetai);
            if and(and(-(1e-18)<=pxij,pxij<=(1+1e-18)),...
                    and(-(1e-18)<=pyij,pyij<=(1+1e-18)))
                p=[p [pxij;pyij]];
                N(i,j+1)=size(p,2);
            end
        end
    end
    p=[p [x2(i);y2(i)]];


    N(i,end)=size(p,2);
end

e=[];
for i=1:n % n = 183
    index1=setdiff(1:(n+2),find(isnan(N(i,:)))); 
    nodes=N(i,index1);
    if abs(p(1,nodes(1))-p(1,nodes(end)))<tol
        %vertical line
        [~,index2]=sort(p(2,nodes));
    else
        [~,index2]=sort(p(1,nodes));
    end

    nodesord=nodes(index2);
    e=[e [nodesord(1:end-1);nodesord(2:end)]];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xstart,xstop,ystart,ystop] = plp(n,tol)
%(poisson line process) but we specify a fix n. Random line coordinates on 
%the boundary of the unit square 

r=1/sqrt(2); %radius of disk
xx0=.5; yy0=.5; %center of disk
%massLine=2*pi*r*lambda; %total measure of the line process

%numbLines=poissrnd(massLine);%Poisson number of lines
numbLines=n; %n lines hard coded
theta=2*pi*rand(numbLines,1); %choose angular component uniformly
p=r*rand(numbLines,1); %choose radial component uniformly
q=sqrt(r.^2-p.^2); %distance to circle edge (alonge line)

%calculate segment endpoints of Poisson line process
xx1=xx0+p.*cos(theta)+q.*sin(theta);
yy1=yy0+p.*sin(theta)-q.*cos(theta);
xx2=xx0+p.*cos(theta)-q.*sin(theta);
yy2=yy0+p.*sin(theta)+q.*cos(theta);

%find start and stop coordinates on unit square
xstart=zeros(size(xx1));xstop=xstart;ystart=xstart;ystop=xstart;
for i=1:length(xx1)
    if or(or(and(xx1(i)<0,xx2(i)<0),and(xx1(i)>1,xx2(i)>1)),...
            or(and(yy1(i)<0,yy2(i)<0),and(yy1(i)>1,yy2(i)>1)))
        %does not cut the unit square 
         xstart(i)=NaN;ystart(i)=NaN;
         xstop(i)=NaN;ystop(i)=NaN;
    else
        %cuts the unit square
        if abs(xx2(i)-xx1(i))<tol
            %vertical line
            xstart(i)=xx1(i);ystart(i)=0;
            xstop(i)=xx1(i);ystop(i)=1;
        else
            %x=x0 when y=0, x=x1 when y=1, y=y0 when x=0, y=y1 when x=1
            theta=atan((yy2(i)-yy1(i))./(xx2(i)-xx1(i)));
            x0=xx1(i)+(0-yy1(i))/tan(theta); 
            y0=yy1(i)+(0-xx1(i))*tan(theta); 
            x1=xx1(i)+(1-yy1(i))/tan(theta);
            y1=yy1(i)+(1-xx1(i))*tan(theta);
            %6 cases
            if and(and(x0>0,x0<1),and(y0>0,y0<1)) %SW
                xstart(i)=x0;ystart(i)=0;
                xstop(i)=0;ystop(i)=y0;
            elseif and(and(x0>0,x0<1),and(x1>0,x1<1)) %SN
                xstart(i)=x0;ystart(i)=0;
                xstop(i)=x1;ystop(i)=1;
            elseif and(and(x0>0,x0<1),and(y1>0,y1<1)) %SE
                xstart(i)=x0;ystart(i)=0;
                xstop(i)=y1;ystop(i)=1;
            elseif and(and(y0>0,y0<1),and(y1>0,y1<1)) %WE
                xstart(i)=0;ystart(i)=y0;
                xstop(i)=1;ystop(i)=y1;
            elseif and(and(x1>0,x1<1),and(y1>0,y1<1)) %NE
                xstart(i)=x1;ystart(i)=1;
                xstop(i)=1;ystop(i)=y1;
            else                                   %NW
                xstart(i)=0;ystart(i)=y0;
                xstop(i)=x1;ystop(i)=1;
            end
        end
    end
end
%assure NaN is not included
indexN=setdiff(1:length(xstart),find(isnan(xstart)));
xstart=xstart(indexN);ystart=ystart(indexN);xstop=xstop(indexN);ystop=ystop(indexN);
%plot([xx1';xx2'],[yy1';yy2'],'r');
end

function [xstart,xstop,ystart,ystop,mids] = fl(nlines,len)
    % Written by Moritz Hauck
    % generate start and stop points of finite length fibers
    xstart = zeros(nlines,1);
    ystart = zeros(nlines,1);
    xstop = zeros(nlines,1);
    ystop = zeros(nlines,1);
    mids = zeros(nlines,2);
    indline = 1;
    while indline <= nlines
        mid = -.5*len + (1+len)*rand(1,2);
        mids(indline,:) = mid;
        phi = pi*rand;
        xmid = mid(1);
        ymid = mid(2);
        xy = [(tan(phi)*xmid-ymid)/tan(phi), 0;...
        (1+tan(phi)*xmid-ymid)/tan(phi), 1;...
        0, ymid-tan(phi)*xmid;...
        1, tan(phi)+ymid-tan(phi)*xmid];
        v = find(all(abs(xy-.5)<=.5,2));
        if isempty(v)
            continue;
        end
        assert(length(v) == 2,'length must be 2.');
        xstart(indline) = xy(v(1),1);
        ystart(indline) = xy(v(1),2);
        xstop(indline) = xy(v(2),1);
        ystop(indline) = xy(v(2),2);
        indline = indline+1;
    end % while
end % function


function [p,e] = discretizefl(xstart,xstop,ystart,ystop,mids,len)
    % Written by Moritz Hauck
    % generate spatial network
    % p matrix
    nlines = length(xstart);
    N=sparse(nlines,nlines+2);p=[];
    for i=1:nlines
        if sum((mids(i,:)-[xstart(i) ystart(i)]).^2) <= len^2
            p=[p [xstart(i);ystart(i)]];
            N(i,1)=size(p,2);
        end % if
        % old nodes
        for j=1:i-1
            N(i,j+1)=N(j,i+1);
        end % for
        for j=i+1:nlines
            % compute angles 
            thetai=atan((ystop(i)-ystart(i))./(xstop(i)-xstart(i)));
            thetaj=atan((ystop(j)-ystart(j))./(xstop(j)-xstart(j)));
            % coord where/if line i and j cut
            R=[cos(thetaj) -cos(thetai);sin(thetaj) -sin(thetai)];
            if abs(det(R))>1e-11 %not parallell
                s=R\[xstart(i)-xstart(j);ystart(i)-ystart(j)];
                pxij=xstart(i)+s(2)*cos(thetai);
                pyij=ystart(i)+s(2)*sin(thetai);
                if sum(([pxij pyij]-mids(i,:)).^2) <= len^2 && sum(([pxij pyij]-mids(j,:)).^2) < len^2 && 0<pxij && pxij<1 && 0<pyij && pyij<1
                    p=[p [pxij;pyij]];
                    N(i,j+1)=size(p,2);
                end % if
            end % if
        end % for
        if sum((mids(i,:)-[xstop(i) ystop(i)]).^2) <= len^2
            p=[p [xstop(i);ystop(i)]];
            N(i,end)=size(p,2);
        end % if
    end % for
    
    % e matrix
    e=[];
    for i=1:nlines
        nodes=nonzeros(N(i,:)); 
        if ~isempty(nodes)
            if abs(p(1,nodes(1))-p(1,nodes(end)))<1e-8
                % vertical line
                [~,indexstop]=sort(p(2,nodes));
            else
                [~,indexstop]=sort(p(1,nodes));
            end % if
            nodesord=nodes(indexstop)';
            e=[e [nodesord(1:end-1);nodesord(2:end)]];
        end % if
    end % for
end % function


