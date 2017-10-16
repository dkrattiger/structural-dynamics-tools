function [i_ss_dofs] = ...
    partition_unit_cell(coordinates,R,K,i_i,i_b,max_ss_nodes,plot_parts)
    
% Dimitri Krattiger (4-3-2017)
%
% description
% ===========
% This code partitions a unit cell for use with the  Automated Multi-Level 
% Substructuring (AMLS) algorithm. The partitioning is performed using the
% spectral algorithm.
%
% inputs
% ======
% coordinates   = node coordinates
% 
% R             = Lattice vectors: [r1,r2,r3]
%
% K             = Stiffness matrix
%
% i_i           = interior DOF index vector
%
% i_b           = boundary DOF index vector
%
% max_ss_nodes  = maximum number of nodes in a partition (substructure)
%
% plot_parts    = true or false, plot partitioned unit cell?
%
% outputs
% =======
% i_ss_dofs     =  cell array tree containing substructure partition
%                  degrees of freedom (for use with AMLS algorithm)

% make sure interior and boundary sets are vectors
i_i = i_i(:);
i_b = i_b(:);

% number of nodes and dimension of problem geometry
[n_nodes,n_dim] = size(coordinates);
n_dof = size(K,1);
n_dpn = n_dof/n_nodes;

% at this stage, need to assume that every node has the same number of DOFs
if rem(size(K,1),n_nodes) ~= 0
    error('no. of deg. of freedom must be divisible by no. of nodes')
end

% connectivity matrix
A = K(1:n_dpn:end,1:n_dpn:end);
for i = 2:n_dpn-1
    A = A + K(i:n_dpn:end,i:n_dpn:end);
end

% perform spectral partitioning (nodes are expressed in lattice vector
% coordinates)

% Not SURE ABOUT THIS
i_i_temp = i_i;


%
i_b_temp = (1:n_nodes)'; i_b_temp(i_i_temp) = [];
[i_ss_nodes_i,i_ss_nodes_b] = spectral_partition(A,max_ss_nodes,i_i_temp,...
    i_b_temp,coordinates/R');

% reorganize substructures so ith boundary is division nodes of ith
% substructure
for i = 1:length(i_ss_nodes_i)
    for j = 1:length(i_ss_nodes_i{i})
        if i<length(i_ss_nodes_i)
            i_ss_nodes_b{i}{j} = i_ss_nodes_i{i}{j};
            i_ss_nodes_b{i}{j} = setdiff(i_ss_nodes_b{i}{j},i_ss_nodes_i{i+1}{(j-1)*2+1});
            i_ss_nodes_b{i}{j} = setdiff(i_ss_nodes_b{i}{j},i_ss_nodes_i{i+1}{(j-1)*2+2});
        else
            i_ss_nodes_b{i}{j} = [];
        end
    end
end

%% plot unit cell paritions
% ======================================================================= %

if plot_parts
    
    figure(2);hold on
    title('Automated Multi Level Substructuring Division')
    
    colors = get(gca,'colororder');         
    markers = {'o','s','*','s','<','>'};
    i_color = 1;
    i = length(i_ss_nodes_i);
    
    if n_dim == 2
        exp_fac = 0.0;
    elseif n_dim == 3
        exp_fac = 0.25;
        exp_fac = 0;
    end
    
    for i = 1:length(i_ss_nodes_i)
        for j = 1:length(i_ss_nodes_i{i})
            
            plot_interior = false;
            if i == length(i_ss_nodes_i)
                plot_interior = true;
            elseif isempty(i_ss_nodes_i{i+1}{2*j}) & isempty(i_ss_nodes_i{i+1}{2*j-1})
                plot_interior = true;
            end            

            if plot_interior
                exp_vec = exp_fac*mean(coordinates(i_ss_nodes_i{i}{j},:));
                c_index = rem(i_color-1,size(colors,1))+1;

                % exploded substructure coordinates and connectivity matrix
                coordinates_ss_exp = coordinates(i_ss_nodes_i{i}{j},:) + ...
                        ones(length(i_ss_nodes_i{i}{j}),1)*exp_vec;
                A_ss = A(i_ss_nodes_i{i}{j},i_ss_nodes_i{i}{j});

                if n_dim == 2
                        
                    if size(coordinates_ss_exp,1)>=3
                        k_bnd = boundary2(coordinates_ss_exp,A_ss)';

                        h1 = patch('faces',k_bnd,'vertices',coordinates_ss_exp);
                        set(h1,'facecolor',colors(c_index,:))
                        set(h1,'edgecolor',colors(c_index,:)*0.85)
                    end

                    plot(coordinates_ss_exp(:,1),coordinates_ss_exp(:,2),...
                         '.','color',colors(c_index,:)*0.85)
                elseif n_dim == 3

                    k_bnd = boundary(coordinates(i_ss_nodes_i{i}{j},1) + exp_vec(1),...
                                      coordinates(i_ss_nodes_i{i}{j},2) + exp_vec(2),...
                                      coordinates(i_ss_nodes_i{i}{j},3) + exp_vec(3));

                    h1 = patch('faces',k_bnd,'vertices',...
                        coordinates(i_ss_nodes_i{i}{j},:) + ...
                            (exp_vec(ones(length(i_ss_nodes_i{i}{j}),1),:)));
                    set(h1,'facecolor',colors(c_index,:))
                    set(h1,'edgecolor',colors(c_index,:)*0.85)

                    plot3(coordinates(i_ss_nodes_i{i}{j},1) + exp_vec(1),...
                          coordinates(i_ss_nodes_i{i}{j},2) + exp_vec(2),...
                          coordinates(i_ss_nodes_i{i}{j},3) + exp_vec(3),...
                          '.','color',colors(c_index,:)*0.85)
                end
                i_color = i_color+1;
            elseif ~plot_interior
                
                coordinates_ss_exp = coordinates(i_ss_nodes_b{i}{j},:)*(1+exp_fac);
                if n_dim == 2
                    
                    xlocs_ss = coordinates_ss_exp(:,1);
                    ylocs_ss = coordinates_ss_exp(:,2);
                    plot(xlocs_ss,ylocs_ss,'k.','markersize',10)


                elseif n_dim == 3

                    b_pts = coordinates(i_ss_nodes_b{i}{j},:);

                    % if x,y, or z-coordinate doesn't change, remove that coordinate
                    keep_col = true(1,3);

                    if size(b_pts,1)>1
                        for k = 1:3 
                            if all(b_pts(:,k) == b_pts(1,k))
                                keep_col(k) = false;
                            end
                        end
                    end
                    b_pts = b_pts(:,keep_col);

                    if ~isempty(b_pts)

                        k_bnd = boundary(b_pts,1);
                        if size(k_bnd,2) == 1;
                            k_bnd = k_bnd';
                        end

                        h1 = patch('faces',k_bnd,'vertices',...
                            coordinates(i_ss_nodes_b{i}{j},:)*(1+exp_fac));
                        set(h1,'facecolor','k')
                        set(h1,'facealpha',0.25)
                        set(h1,'edgecolor','k')

                    end
                end
            end
        end
    end
            
    axis equal
    if n_dim == 3
        view(3)
    end
    drawnow
end

%% create single index for substructures
n_lev = length(i_ss_nodes_i);
for i = 1:n_lev
    for j = 1:length(i_ss_nodes_i{i})
        % dofs in substructure i
        if i == n_lev
            i_ss_nodes{i}{j} = i_ss_nodes_i{i}{j};
        else
            i_ss_nodes{i}{j} = i_ss_nodes_b{i}{j};
        end
    end
end

% shift all substructure down one level in substructure tree
n_lev = length(i_ss_nodes);
for i = 1:n_lev
    i_ss_nodes_new{i+1} = i_ss_nodes{i};
end

% add boundary dofs as root level substructure
i_ss_nodes_new{1}{1} = i_b;

% convert node list to DOF list
for i = 1:length(i_ss_nodes_new)
    for j = 1:length(i_ss_nodes_new{i})
        i_ss_dofs{i}{j} = node2dof(i_ss_nodes_new{i}{j},n_dpn);
    end
end