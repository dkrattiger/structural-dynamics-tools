function [K_AMLS,M_AMLS,n_FI_ss,varargout] = AMLS(i_ss_dofs,K_free,M_free,w_cut,dof_sets,varargin)
%
% Dimitri Krattiger
%
% Description
% ===========
% This code computes the Automated Multi Level Substructure (AMLS) 
% representation mass and stiffness matrices. This is essentially a 
% recursive Craig-Bampton (CB) calculation that should be faster than a 
% single-level CB calculation. The main difference from CB to AMLS is that
% interfaces between substructures are not treated again at every level,
% rather they are treated like they are substructures themselves.
%
% Inputs
% ======
% i_ss_dofs = tree shaped cell array containing DOF indices of partitioned
%             substructures
%
% K_free    = free stiffness matrix
% 
% M_free    = free mass matrix
%
% w_cut     = cutoff frequency for fixed-interface mode calculations
% 
% options   = structure containing optional parameters
%
%
% Outputs
% =======
% K_AMLS    = AMLS stiffness matrix
% 
% M_AMLS    = AMLS mass matrix
%
% n_FI_ss   = number of fixed-interface modes from each substructure
%
% T_AMLS    = transformation between full DOF vector and AMLS DOF 
%             vector
% 
% Frs       = residual flexibility matrix
%             (only generated if options.resPlus = true)
%
% Psi_hat   = constraint mode matrix
%             (only generated if options.resPlus = true)

%% Check what inputs are given and set rest to default
% ======================================================================= %

% specify default options
defaults.verbose        = false;
defaults.verboseTab     = '%%   ';
defaults.verboseTabSum  = '';
defaults.resPlus        = false;
defaults.outputT        = false;

if nargin>4
    options = varargin{1};
    options = setstructfields(defaults,options);
else
    options = defaults;
end

%% Display Timing info
% ======================================================================= %

if options.verbose
    t_start_AMLS = tic;    
    tb = options.verboseTabSum;
    options.verboseTabSum = [options.verboseTabSum,options.verboseTab];
    tbi = options.verboseTabSum;
    fprintf([tb,'\n'])
    fprintf([tb,'Automated Multi Level Substructuring (AMLS)\n'])
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
end

%% create tree-to-number cell array
count = 0;
% for i = 1:length(i_ss_dofs)
for i = length(i_ss_dofs):-1:1
    for j = 1:length(i_ss_dofs{i})
        count = count + 1;
        tree2num{i}{j} = count;
        num2tree(count,:) = [i,j];        
    end
end
n_ss = count;

%% create ancestor tree
% ======================================================================= %
ancestors{1}{1} = [0,1];
ancestors{1}{1} = [];
ancestors{2}{1} = [1,1];
for i = 3:length(i_ss_dofs)
    for j = 1:length(i_ss_dofs{i})
        ancestors{i}{j} = [ancestors{i-1}{ceil(j/2)};[i-1,ceil(j/2)]];
    end
end

%% create descendent tree;
% ======================================================================= %
descendents{1}{1} = [0,1];
descendents{1}{1} = [];
descendents{2}{1} = [1,1];
for i = length(i_ss_dofs):-1:2
    for j = 1:length(i_ss_dofs{i})
        if i == length(i_ss_dofs)
            descendents{i}{j} = [];
        else
            descendents{i}{j} = [i+1,2*(j-1)+1;descendents{i+1}{2*(j-1)+1};...
                                 i+1,2*(j-1)+2;descendents{i+1}{2*(j-1)+2}];
        end
    end
end
descendents{1}{1} = [[2,1];descendents{2}{1}];

%% split mass and stiffness into cell array
% ======================================================================= %
% this will allow us to change the matrix block by block rather than all at
% once
K_CB_cell = cell(n_ss,n_ss);
M_CB_cell = cell(n_ss,n_ss);

% add diagonal entries into mass and stiffness for all but bottom level
% (from the bottom level these appear as boundary-boundary entries

for i = 1:length(i_ss_dofs)
    for j = 1:length(i_ss_dofs{i})
        substr_ind = tree2num{i}{j};
    
        % dofs in substructure i
        i_ii = i_ss_dofs{i}{j};
        n_dof_ss(substr_ind) = length(i_ii);
        
        
    
        % number of ancestors for current substructure
        n_ancestors = size(ancestors{i}{j},1);

        for k = 1:n_ancestors
            
            ancestor_lev = ancestors{i}{j}(k,1);
            ancestor_num = ancestors{i}{j}(k,2);
            ancestor_ind = tree2num{ancestor_lev}{ancestor_num};
            
            % dofs in ancestor k
            i_kk = i_ss_dofs{ancestor_lev}{ancestor_num};
            
            % off diagonal mass and stiffness blocks
            K_CB_cell{substr_ind,ancestor_ind} = K_free(i_ii,i_kk);
            M_CB_cell{substr_ind,ancestor_ind} = M_free(i_ii,i_kk);
        end
        
        % diagonal mass and stiffness blocks
        K_CB_cell{substr_ind,substr_ind} = K_free(i_ii,i_ii);
        M_CB_cell{substr_ind,substr_ind} = M_free(i_ii,i_ii);

    end
end

%% perform mass and stiffness reductions
% ======================================================================= %
% preallocate constraint mode array
Psi_cell = cell(n_ss,n_ss);
PHI_cell = cell(n_ss,1);
if options.resPlus
    Frs_cell = cell(n_ss,1);
end

for i = length(i_ss_dofs):-1:2
    t_start_level = tic;
    
    % number of substructures in current tree level
    n_leafs = length(i_ss_dofs{i});
    
    % display run info
    if options.verbose
        fprintf([tbi,'LEVEL %i, %i LEAFS\n'],i,n_leafs)
        fprintf([tbi,'==========================\n'])
    end
        
    for j = 1:n_leafs
        t_start_leaf = tic;
        
        % substructure indices (tree notation and index notation)
        substr_ind = tree2num{i}{j};
        
        
        full_eig = n_dof_ss(substr_ind) <= 800;
        if full_eig
            
            [PHI_FI,L] = eig(full(K_CB_cell{substr_ind,substr_ind}),...
                full(M_CB_cell{substr_ind,substr_ind}));

        else 
            % use iterative eigenvalue solver to compute fixed interface
            % modes if substructure is large
            
            % initialize eigenvalue trial loop
            n_eigs = 10;
            L = 0;
            
            % increase number of modes computed until it spans the proper
            % range
            while max(max(L))<w_cut^2 && n_eigs < size(K_CB_cell{substr_ind,substr_ind},1)
                    
                % compute eigenvalues
                [PHI_FI,L] = eigs(K_CB_cell{substr_ind,substr_ind},...
                    M_CB_cell{substr_ind,substr_ind},n_eigs,'sm');                
                
                % double number of eigenvalues to compute next time through
                n_eigs = 2*n_eigs;
            end
            
            % check if we should have done full eig solution after all
            if max(max(L))<w_cut^2 && n_eigs > size(K_CB_cell{substr_ind,substr_ind},1)
                full_eig = true; 
                [PHI_FI,L] = eig(full(K_CB_cell{substr_ind,substr_ind}),...
                            full(M_CB_cell{substr_ind,substr_ind}));
            end
        end
        
        % sort eigenvalues and eigenvectors
        [L,i_sort] = sort(diag(L));
        PHI_FI = PHI_FI(:,i_sort);
        
        % for degenerate modes, perform a secondary diagonalization to
        % ensure that modes are mass orthogonal
        L_diff = L(2:end)-L(1:end-1);
        if any(L_diff<1e-12*w_cut)
            i_deg = find(L_diff<1e-12*w_cut);
            i_deg = unique([i_deg;i_deg+1]);
            [V,~] = eig(PHI_FI(:,i_deg)'*M_CB_cell{substr_ind,substr_ind}*PHI_FI(:,i_deg));
            PHI_FI(:,i_deg) = PHI_FI(:,i_deg)*V;
        end
        
            
        % mass normalize fixed interface mode shapes
        PHI_FI = PHI_FI*diag(diag(PHI_FI'*M_CB_cell{substr_ind,substr_ind}*PHI_FI).^(-0.5));
        
        % determine number of fixed interface modes
        n_FI = sum(L<=w_cut^2);
        
        % dominant fixed interface modes
        PHId = PHI_FI(:,1:n_FI);
        Ld = L(1:n_FI);
        
        if full_eig
            PHIr = PHI_FI(:,n_FI+1:end);
            Lr = L(n_FI+1:end);
        end
        
        % Fixed interface mode storage array
        PHI_cell{substr_ind} = PHId;        

        % number of ancestors for current substructure
        n_ancestors = size(ancestors{i}{j},1);

        for k = 1:n_ancestors

            % dofs in k-th ancestor
            ancestor_lev = ancestors{i}{j}(k,1);
            ancestor_num = ancestors{i}{j}(k,2);
            ancestor_ind1 = tree2num{ancestor_lev}{ancestor_num};

            % constraint modes
            if isempty(Psi_cell{substr_ind,ancestor_ind1})
                Psi_cell{substr_ind,ancestor_ind1} = ...
                    full(-mldivide(K_CB_cell{substr_ind,substr_ind},...
                                   K_CB_cell{substr_ind,ancestor_ind1}));
            end
            
            Mc = Psi_cell{substr_ind,ancestor_ind1}'*M_CB_cell{substr_ind,substr_ind} + ...
                    M_CB_cell{substr_ind,ancestor_ind1}'; 
            
            Mctest = Psi_cell{substr_ind,ancestor_ind1}'*M_CB_cell{substr_ind,substr_ind} + ...
                    M_CB_cell{substr_ind,ancestor_ind1}'; 
            for p = 1:(k-1)

                % dofs in p-th ancestor
                ancestor_lev = ancestors{i}{j}(p,1);
                ancestor_num = ancestors{i}{j}(p,2);
                ancestor_ind2 = tree2num{ancestor_lev}{ancestor_num};

                if isempty(Psi_cell{substr_ind,ancestor_ind2})
                    Psi_cell{substr_ind,ancestor_ind2} = -mldivide(K_CB_cell{substr_ind,substr_ind},K_CB_cell{substr_ind,ancestor_ind2});
                end

                % compute modification to boundary-boundary mass and stiffness
                m_add = Mc*Psi_cell{substr_ind,ancestor_ind2} + Psi_cell{substr_ind,ancestor_ind1}'*M_CB_cell{substr_ind,ancestor_ind2};
                M_CB_cell{ancestor_ind1,ancestor_ind2} = M_CB_cell{ancestor_ind1,ancestor_ind2} + m_add;

                k_add = K_CB_cell{substr_ind,ancestor_ind1}'*Psi_cell{substr_ind,ancestor_ind2};
                K_CB_cell{ancestor_ind1,ancestor_ind2} = K_CB_cell{ancestor_ind1,ancestor_ind2} + k_add;
            end

            % compute modification to boundary-boundary mass and stiffness
            m_add = Mc*Psi_cell{substr_ind,ancestor_ind1} + Psi_cell{substr_ind,ancestor_ind1}'*M_CB_cell{substr_ind,ancestor_ind1};
            M_CB_cell{ancestor_ind1,ancestor_ind1} = M_CB_cell{ancestor_ind1,ancestor_ind1} + m_add;

            k_add = K_CB_cell{substr_ind,ancestor_ind1}'*Psi_cell{substr_ind,ancestor_ind1};
            K_CB_cell{ancestor_ind1,ancestor_ind1} = K_CB_cell{ancestor_ind1,ancestor_ind1} + k_add;
            modified(ancestor_ind1,ancestor_ind1) = 1;
        end

        % compute enhanced modification to boundary-boundary mass and
        % stiffness
        if options.resPlus
            if full_eig
                Frs_cell{substr_ind} = PHIr*diag(Lr.^(-1))*PHIr';
            else
                Frs_cell{substr_ind} = inv(K_CB_cell{substr_ind,substr_ind})-PHId*diag(Ld.^(-1))*PHId';
            end        
        end
        
        % use fixed interface modes to transform off-diagonal terms in mass and
        % stiffness matrices
        for k = 1:n_ancestors

            % dofs in k-th ancestor
            ancestor_lev = ancestors{i}{j}(k,1);
            ancestor_num = ancestors{i}{j}(k,2);
            ancestor_ind1 = tree2num{ancestor_lev}{ancestor_num};

            % off-diagonal mass term (off diagonal stiffness term goes to 0)
            mu = PHId'*(M_CB_cell{substr_ind,substr_ind}*Psi_cell{substr_ind,ancestor_ind1}+...
                M_CB_cell{substr_ind,ancestor_ind1});
            M_CB_cell{substr_ind,ancestor_ind1} =  mu;
            K_CB_cell{substr_ind,ancestor_ind1} =  zeros(size(mu));
            
        end

        % add diagonal terms into mass and stiffness matrix
        K_CB_cell{substr_ind,substr_ind} = diag(Ld);
        M_CB_cell{substr_ind,substr_ind} = eye(n_FI);
        modified(substr_ind,substr_ind) = 1;
        
        % number of descendents for current substructure
        n_descendents = size(descendents{i}{j},1);
        
        % add off-diagonal descendent terms
        for k = 1:n_descendents
            
            % dofs in k-th ancestor
            descendent_lev = descendents{i}{j}(k,1);
            descendent_num = descendents{i}{j}(k,2);
            descendent_ind1 = tree2num{descendent_lev}{descendent_num};
            
            for p = 1:n_ancestors
                
                % dofs in p-th ancestor
                ancestor_lev = ancestors{i}{j}(p,1);
                ancestor_num = ancestors{i}{j}(p,2);
                ancestor_ind1 = tree2num{ancestor_lev}{ancestor_num};
                                
                mu = M_CB_cell{descendent_ind1,ancestor_ind1} + ...
                     M_CB_cell{descendent_ind1,substr_ind}*Psi_cell{substr_ind,ancestor_ind1};

                
                M_CB_cell{descendent_ind1,ancestor_ind1} = mu;
                K_CB_cell{descendent_ind1,ancestor_ind1} = zeros(size(mu));
            end            
            
            MA = M_CB_cell{descendent_ind1,substr_ind}*PHId;
            M_CB_cell{descendent_ind1,substr_ind} = MA;
            K_CB_cell{descendent_ind1,substr_ind} = zeros(size(MA));
        end
        
        % display timing info
        if options.verbose
            if full_eig
                fprintf([tbi,'\tLeaf %i of %i:\t%iDOFs,\t%6.2fsec (direct eigenvalue solution)\n'],...
                    j,n_leafs,n_dof_ss(substr_ind),toc(t_start_leaf))
            else
                fprintf([tbi,'\tLeaf %i of %i:\t%iDOFs,\t%6.2fsec (iterative eigenvalue solution)\n'],...
                    j,n_leafs,n_dof_ss(substr_ind),toc(t_start_leaf))
            end
        end
            
    end
    
    if options.verbose
        fprintf([tbi,'\tall %i leafs: %6.2fsec \n'],n_leafs,toc(t_start_level))
    end
end

%% Direct Mass and stiffness assembly
% ======================================================================= %
tic
for i = 1:n_ss
    n_FI_ss(i) = size(M_CB_cell{i,i},1);
end

% initialize sparse matrix formation vectors
n_dof_AMLS = sum(n_FI_ss);
K_AMLS = zeros(n_dof_AMLS,n_dof_AMLS);
M_AMLS = zeros(n_dof_AMLS,n_dof_AMLS);

for i = 1:length(i_ss_dofs)
    for j = 1:length(i_ss_dofs{i})
        substr_ind = tree2num{i}{j};
        n_row_start = sum(n_FI_ss(1:(substr_ind-1)));
        i_row = (n_row_start+1):n_row_start+n_FI_ss(substr_ind);
        
        % diagonal stiffness and mass blocks
        K_AMLS(i_row,i_row) = K_CB_cell{substr_ind,substr_ind};
        M_AMLS(i_row,i_row) = M_CB_cell{substr_ind,substr_ind};
        
        % number of ancestors for current substructure
        n_ancestors = size(ancestors{i}{j},1);

        for k = 1:n_ancestors
            
            % index of ancestor
            ancestor_lev = ancestors{i}{j}(k,1);
            ancestor_num = ancestors{i}{j}(k,2);
            ancestor_ind = tree2num{ancestor_lev}{ancestor_num};
            
            % start of column index
            n_col_start = sum(n_FI_ss(1:(ancestor_ind-1)));
            i_col = (n_col_start+1):n_col_start+n_FI_ss(ancestor_ind);
            
            % off-diagonal mass blocks
            M_AMLS(i_row,i_col) = M_CB_cell{substr_ind,ancestor_ind};
            M_AMLS(i_col,i_row) = M_CB_cell{substr_ind,ancestor_ind}';

        end
    end
end

if options.verbose
    fprintf([tbi,'Mass, Stiffness direct assembly time: %6.2f \n'],toc)
end

%% Form AMLS transformation directly
% ======================================================================= %

if options.outputT | options.resPlus
    
% %% Form AMLS transformation directly
% ======================================================================= %
% % number of fixed interface modes in each substructure
% n_FI_ss = zeros(1,n_ss);
% for i = 1:n_ss
%     n_FI_ss(i) = size(M_CB_cell{i,i},1);
% end

    % last entry of PHI_cell is for the top level boundary. Make this
    % identity because boundary reduction will occur later
    PHI_cell{end} = eye(n_FI_ss(end));  % n_dof_ss(end)  should=  n_FI_ss(end)
    if options.resPlus
        Frs_cell{end} = zeros(n_FI_ss(end),n_FI_ss(end));
    end
    
    
    % number of dofs in original model
    n_dof = sum(n_dof_ss);

    % % number of dofs in AMLS reduced model
    % n_dof_AMLS = sum(n_FI_ss);

    % allocate space for arrays used in forming enhanced AMLS transformation
    T_AMLS = zeros(n_dof,n_dof_AMLS);
    Frs = spalloc(n_dof,n_dof,sum(n_dof_ss.^2));

    n_el_Psi = 0;
    for i = 1:size(Psi_cell,1)
        for j = 1:size(Psi_cell,2)
            n_el_Psi = n_el_Psi + numel(Psi_cell{i,j});
        end
    end
    n_el_Psi = n_el_Psi + n_dof;

    % preallocate Psi_hat array sparse triplet vectors
    Psi_hat_rows = zeros(n_el_Psi,1);
    Psi_hat_cols = zeros(n_el_Psi,1);
    Psi_hat_vals = zeros(n_el_Psi,1);

    % initialize counter
    count = 0;

    % Assemble global constraint mode matrix
    % NOTE THAT Psi_cell IS DESTROYED BY THIS PROCESS
    for i = length(i_ss_dofs):-1:1
        for j = 1:length(i_ss_dofs{i})
            substr_ind = tree2num{i}{j};

            % dofs in substructure i (original sort)
            i_rows = i_ss_dofs{i}{j};

            % number of ancestors for current substructure
            n_ancestors = size(ancestors{i}{j},1);

            for k = n_ancestors:-1:1

                ancestor_lev = ancestors{i}{j}(k,1);
                ancestor_num = ancestors{i}{j}(k,2);
                ancestor_ind = tree2num{ancestor_lev}{ancestor_num};

                % number of ancestors of substructure i that are all
                % descendents of ancestor k
                middle_men = find(ancestors{i}{j}(:,1)>ancestor_lev);
                n_middle_men =  length(middle_men);

                Psi_hat_sub = Psi_cell{substr_ind,ancestor_ind};


                for p = 1:n_middle_men

                    middle_lev = ancestors{i}{j}(middle_men(p),1);
                    middle_num = ancestors{i}{j}(middle_men(p),2);
                    middle_ind = tree2num{middle_lev}{middle_num};

                    Psi_hat_sub = Psi_hat_sub + full(Psi_cell{substr_ind,middle_ind})*full(Psi_cell{middle_ind,ancestor_ind});

                end

                % update Psi cell array
                Psi_cell{substr_ind,ancestor_ind} = Psi_hat_sub;
                i_cols = (sum(n_FI_ss(1:(ancestor_ind-1)))+(1:n_FI_ss(ancestor_ind)))';

                % AMLS transformation matrix
                T_AMLS(i_rows,i_cols) = Psi_hat_sub*PHI_cell{ancestor_ind};

                % index of global columns in which to drop subsstructure
                % constraint modes
                i_cols2 = (sum(n_dof_ss(1:(ancestor_ind-1)))+(1:n_dof_ss(ancestor_ind)))';

                % assemble global constraint mode matrix
                [i_col_grid,i_row_grid] = meshgrid(i_cols2,i_rows);
                n_el = length(i_row_grid(:));            
                Psi_hat_rows(count + (1:n_el)) = i_row_grid(:);
                Psi_hat_cols(count + (1:n_el)) = i_col_grid(:);
                Psi_hat_vals(count + (1:n_el)) = Psi_hat_sub(:);            
                count = count + n_el;
            end

            % dofs in AMLS substructure
            i_cols = (sum(n_FI_ss(1:(substr_ind-1)))+(1:n_FI_ss(substr_ind)))';

            % diagonal blocks of AMLS transformation
            T_AMLS(i_rows,i_cols) = PHI_cell{substr_ind};

            % dofs in substructure (block sort)
            i_cols2 = (sum(n_dof_ss(1:(substr_ind-1)))+(1:n_dof_ss(substr_ind)))';

            % global constraint mode matrix
            n_el = n_dof_ss(substr_ind);
            Psi_hat_rows(count + (1:n_el)) = i_rows(:);
            Psi_hat_cols(count + (1:n_el)) = i_cols2(:);
            Psi_hat_vals(count + (1:n_el)) = ones(n_dof_ss(substr_ind),1);
            count = count + n_el;

            % add substructure residual flexibility into global residual
            % flexibility matrix
            if options.resPlus
                Frs(i_cols2,i_cols2) = Frs_cell{substr_ind};
            end

        end
    end
    
    if options.verbose
        fprintf([tbi,'T_AMLS assembly time: %6.2f \n'],toc)
    end
    
else
    T_AMLS =[];
end
% NOTE THAT Psi_cell IS DESTROYED BY THIS PROCESS

varargout{1} = T_AMLS;

%% Residual Enhanced outputs
if options.resPlus
    
    varargout{2} = Frs;
    
    % form sparse constraint mode matrix from sparse triplet vectors
    tic
    Psi_hat = sparse(Psi_hat_rows,Psi_hat_cols,Psi_hat_vals,n_dof,n_dof);
    if options.verbose
        fprintf([tbi,'Psi_hat sparse matrix formation call: t = %4.2f\n'],toc);
    end
   
    varargout{3} = Psi_hat;
end

%% Find updated dof_sets structure
% ======================================================================= %

dof_set_fields = fieldnames(dof_sets);
dof_sets_mod.i = (1:sum(n_dof_ss(1:end-1)))';
count = max(dof_sets_mod.i);

% boundary DOFs are sorted by boundary set
for i = 1:length(dof_set_fields{i})
    if ~isequal(dof_set_fields{i},'i');
        n_dof_set = length(dof_sets.(dof_set_fields{i}));
        dof_sets_mod.(dof_set_fields{i}) = count + (1:n_dof_set)';
        count = count + n_dof_set;
    end
end

%% display timing info
% ======================================================================= %
if options.verbose
    fprintf([tbi,'AMLS calculation time: %5.2f s\n'],toc(t_start_AMLS))
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
    fprintf([tb,'\n'])
end