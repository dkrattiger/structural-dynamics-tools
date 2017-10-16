function [i_ss_i,i_ss_b] = spectral_partition(A,max_ss_nodes,...
    i_i,i_b,coordinates)
    
% Dimitri Krattiger (4-3-2017)
%
% description
% ===========
% This code uses the spectral algorithm to bisect the
% node sets enough times such that substructure size is less than
% "max_ss_nodes".
%
% inputs
% ======
% A             = Adjacency matrix (specifying which nodes are linked.
%
% max_ss_nodes  = maximum substructure size
%
% i_i           = interior node index
%
% i_b           = boundary node index
%
% max_ss_nodes  = maximum number of nodes in a partition (substructure)
%
% outputs
% =======
% i_ss_i        = tree structure containing interior node indices of every
%                 substructure
%
% i_ss_b        = tree structure containing boundary node indices of every
%                 substructure


% number of nodes and number of dimensions
[n_nodes,n_dim] = size(coordinates);

% show diagnostic plots?
visualize = false;

% store base level info into cell array
i_i = i_i(:);
i_b = i_b(:);
i_ss_i{1}{1} = i_i;
i_ss_b{1}{1} = i_b;
i_ss{1}{1} = sort([i_i;i_b]);

% initialize substrucure size tracker
ss_size = length(i_ss_i{1}{1});
i = 1;
while ss_size > max_ss_nodes
    
    % initialize domain plot figure
    if visualize
        hfig = figure(500+i);clf
        
        % use this to plot over the top of a unit cell diagram
        if false
            copyobj(get(2,'children'),hfig); 
        end
        
        colors = get(gca,'colororder');
        markers = {'o','s','*','s','<','>'};
        i_color = 1;
    end
    
    ss_size = 0;
    for j = 1:2^(i-1)
        
        if length(i_ss_i{i}{j}) < max_ss_nodes
            i_ss_i{i+1}{(j-1)*2+1} = [];
            i_ss_b{i+1}{(j-1)*2+1} = [];
            i_ss{i+1}{(j-1)*2+1} = [];

            i_ss_i{i+1}{(j-1)*2+2} = [];
            i_ss_b{i+1}{(j-1)*2+2} = [];
            i_ss{i+1}{(j-1)*2+2} = [];
        else
            
            % pick up interior and boundary partitions from base segment (will
            % be used to identify interior and boundary partitions of split
            % segments)
            i_i = i_ss_i{i}{j};
            i_b = i_ss_b{i}{j};

            % Form Adjacency matrix for current segment
            G = A(i_ss{i}{j},i_ss{i}{j});

            % form laplacian matrix from stiffness 
            L = - spones(G|G');
            L = L - diag(sum(L));        
            L = (1/2)*(L+L');

            % Typically the spectral algorithm uses only the 2nd-lowest
            % eigenvector (the Fiedler vector) to perform partitioning. Due
            % to domain symmetries however, the Fiedler vector may be
            % degenerate and the resulting partition will be due to a
            % random linear combination of the degenerate vectors. Since we
            % would like replicable results, this code tries to form a
            % linear combination of the degenerate vectors that aligns with
            % the coordinate axes.
            
            % compute lowest six eigenpairs
            [V,D] = eigs(L,6,'sm');
            [D,i_sort] = sort(diag(D));
            V = V(:,i_sort); 
            
            % if Fiedler vector is degenerate, use a combination 
            % of degenerate vectors that is parallel to a lattice vector
            degen_tol = 1e-6;
            i_keep = find(abs(D-D(2)) < degen_tol*max(D));
            V = V(:,i_keep);
            [V,~] = qr(V,0);
            range = max(coordinates(i_ss{i}{j},:))-min(coordinates(i_ss{i}{j},:));
            [~,i_sort] = sort(range,'descend');
            if size(V,2) >= 1
                rhsvec = zeros(size(V,2),1);
                rhsvec(1) = 1;
                cvec = (coordinates(i_ss{i}{j},i_sort(1:size(V,2)))'*V)\rhsvec;
                V = V*cvec;
            else
                V = sum(V,2);
            end 
            [Vsort,i_sort] = sort(V);

            % Choose which nodes are in each part by finding jumps in
            % eigenvector values that split eigenvector approximately in
            % half. The weight is somewhat arbitrary and is designed to favor
            % splittings that divide the nodes more evenly in half.
            jumps = (Vsort(2:end)-Vsort(1:end-1))';
            jumps = jumps/max(jumps);
            weight = (1-(([1:length(V)-1]-length(V)/2).^2)/(length(V)/2)^2).^(10);
            split_factor = jumps.*weight;
            [~,i_split] = max(split_factor); 


            % re-sort node index to correspond to original sorting
            if i_split >= length(V)/2
                part1 = (i_sort((i_split+1):length(V)));
                part2 = (i_sort(1:i_split));
            else
                part1 = (i_sort(1:i_split));
                part2 = (i_sort((i_split+1):length(V)));
            end

            % visualize splitting heuristic
            if visualize
                figure(200);clf
                h1 = plot(weight,'r:');hold on
                h2 = plot(jumps,'b:');
                h3 = plot(split_factor,'k');
                h4 = plot(i_split,split_factor(i_split),'mo');
                legend([h1,h2,h3,h4],'weight','jumps','combined','maximum')
            end

            % map bisection information to global node numbers
            part1 = i_ss{i}{j}(part1);
            part2 = i_ss{i}{j}(part2);   

            % add nodes in part 2 to boundary if they share connectivity with
            % nodes in part 1
            Amod = logical(A(part1,:));
            i_part2 = false(length(part1),n_nodes);
            i_part2(:,part2) = true(length(part1),length(part2));        
            i_split = Amod & i_part2;
            split1 = part1(any(i_split,2));
            split2 = find(any(i_split))';
            split = [split1;split2];
            split = unique(split);

            % remove split nodes from each side
            part1 = setdiff(part1,split);
            part2 = setdiff(part2,split);

            % if boundary nodes share no connectivity with part 2
            % remove them from boundary and add to part 1
            connect = false(size(split));
            for j2 = 1:length(split)
                j3 = split(j2);

                connect(j2) = any(A(j3,part2));
            end
            part1 = [part1;split(~connect)];
            split = split(connect);

            % store bisection information in interior and boundary arrays
            i_ss_i{i+1}{(j-1)*2+1} = intersect(i_i,part1,'stable');
            i_ss_b{i+1}{(j-1)*2+1} = [split;intersect(i_b,part1,'stable')];
            i_ss{i+1}{(j-1)*2+1} = sort([i_ss_i{i+1}{(j-1)*2+1};i_ss_b{i+1}{(j-1)*2+1}]);

            i_ss_i{i+1}{(j-1)*2+2} = intersect(i_i,part2,'stable');
            i_ss_b{i+1}{(j-1)*2+2} = [split;intersect(i_b,part2,'stable')];
            i_ss{i+1}{(j-1)*2+2} = sort([i_ss_i{i+1}{(j-1)*2+2};i_ss_b{i+1}{(j-1)*2+2}]);

            % plot partitioned nodes
            if visualize
                figure(hfig);hold on
                if n_dim == 2

                    plot(coordinates(i_ss_i{i+1}{(j-1)*2+1},1),...
                         coordinates(i_ss_i{i+1}{(j-1)*2+1},2),...
                         '.','color',colors(rem(i_color-1,7)+1,:))
                    plot(coordinates(i_ss_b{i+1}{(j-1)*2+1},1),...
                         coordinates(i_ss_b{i+1}{(j-1)*2+1},2),...
                         'linestyle','none',...
                         'marker',markers{rem(i_color-1,6)+1},...
                         'color',colors(rem(i_color-1,7)+1,:),'linewidth',1.5)
                    i_color = i_color+1;

                    plot(coordinates(i_ss_i{i+1}{(j-1)*2+2},1),...
                        coordinates(i_ss_i{i+1}{(j-1)*2+2},2),...
                        '.','color',colors(rem(i_color-1,7)+1,:))
                    plot(coordinates(i_ss_b{i+1}{(j-1)*2+2},1),...
                        coordinates(i_ss_b{i+1}{(j-1)*2+2},2),...
                        'linestyle','none',...
                        'marker',markers{rem(i_color-1,6)+1},...
                        'color',colors(rem(i_color-1,7)+1,:),'linewidth',1.5)
                    i_color = i_color+1;
                elseif n_dim == 3
                    plot3(coordinates(i_ss_i{i+1}{(j-1)*2+1},1),...
                         coordinates(i_ss_i{i+1}{(j-1)*2+1},2),...
                         coordinates(i_ss_i{i+1}{(j-1)*2+1},3),...
                         '.','color',colors(rem(i_color-1,7)+1,:))
                    plot3(coordinates(i_ss_b{i+1}{(j-1)*2+1},1),...
                         coordinates(i_ss_b{i+1}{(j-1)*2+1},2),...
                         coordinates(i_ss_b{i+1}{(j-1)*2+1},3),...
                         'linestyle','none',...
                         'marker',markers{rem(i_color-1,6)+1},...
                         'color',colors(rem(i_color-1,7)+1,:),'linewidth',1.5)
                    i_color = i_color+1;

                    plot3(coordinates(i_ss_i{i+1}{(j-1)*2+2},1),...
                        coordinates(i_ss_i{i+1}{(j-1)*2+2},2),...
                        coordinates(i_ss_i{i+1}{(j-1)*2+2},3),...
                        '.','color',colors(rem(i_color-1,7)+1,:))
                    plot3(coordinates(i_ss_b{i+1}{(j-1)*2+2},1),...
                        coordinates(i_ss_b{i+1}{(j-1)*2+2},2),...
                        coordinates(i_ss_b{i+1}{(j-1)*2+2},3),...
                        'linestyle','none',...
                        'marker',markers{rem(i_color-1,6)+1},...
                        'color',colors(rem(i_color-1,7)+1,:),'linewidth',1.5)
                    i_color = i_color+1;
                    view(3)
                end

                axis equal
                drawnow
                pause
            end      
            ss_size = max([ss_size,length(i_ss_i{i+1}{(j-1)*2+2}),length(i_ss_i{i+1}{(j-1)*2+1})]);

        end
    end
    
    % increment bisection count index i
    i = i + 1;
end