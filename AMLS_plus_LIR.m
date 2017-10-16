function [K_AMLS,M_AMLS,dof_sets_mod,varargout] = AMLS_plus_LIR(i_ss_dofs,dof_sets,K_free,M_free,w_i,options)
%
% Dimitri Krattiger
%
% Description
% ===========
% This code computes the Automated Multi Level Substructure (AMLS) 
% representation mass and stiffness matrices, and the local interface
% reduction (LIR), and then applies residual mode enhancement to the
% interior reduction (plus).
%
% Inputs
% ======
% i_ss_dofs = tree shaped cell array containing DOF indices of partitioned
%             substructures
%
% dof_sets  = DOF set structure for full model
%
% K_free    = free stiffness matrix
% 
% M_free    = free mass matrix
%
% w_i       = cutoff frequency for fixed-interface mode calculations
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
% dof_sets_mod   = DOF set structure for reduced order model
%
% T_AMLS    = transformation between full DOF vector and AMLS DOF 
%             vector

%% Get options arguments and set defaults
% ======================================================================= %

% set_default values for options
defaults.BoundaryMethod = 'none';
defaults.n_CC           = [];
defaults.w_b            = [];
defaults.verbose        = false;
defaults.verboseTab     = '%%   ';
defaults.verboseTabSum  = '';
defaults.outputT        = false;
defaults.LIRorthotype   = 'qr';

options = setstructfields(defaults,options);

%% Display Timing info
% ======================================================================= %

if options.verbose
    t_start_AMLSpl = tic;    
    tb = options.verboseTabSum;
    options.verboseTabSum = [options.verboseTabSum,options.verboseTab];
    tbi = options.verboseTabSum;
    fprintf([tb,'\n'])
    fprintf([tb,'Residual Enhanced Multi-Level Substructuring (AMLS+)\n'])
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
end

%% perform basic AMLS transform
% ======================================================================= %
options_AMLS         = options;
options_AMLS.outputT = true;
options_AMLS.resPlus = true;
[K_AMLS.w0,M_AMLS.w0,n_FI_ss,T_AMLS.w0,Frs,Psi_hat] = AMLS(i_ss_dofs,K_free,M_free,w_i,dof_sets,options_AMLS);

%% Update dof_sets structure
% ======================================================================= %
% update number of fixed interface modes
n_FI = sum(n_FI_ss(1:end-1));

% Create new dof-sets structure that contains info about interior dofs
% and boundary dof sets

% new degree of freedom sets
dof_set_names = fieldnames(dof_sets);
dof_sets_mod.i = (1:n_FI)';
count = n_FI;

% boundary order hasn't changed so just assign new number sets to the
% dof_sets structure
for i = 2:27
    n_dof_set = length(dof_sets.(dof_set_names{i}));
    dof_sets_mod.(dof_set_names{i}) = count + (1:n_dof_set)';
    count = count + n_dof_set;
end

%% Perform LIR reduction
% ======================================================================= %
if ~strcmpi(options.BoundaryMethod,'none')
        
    % perform LIR reduction
    [K_AMLS.w0,M_AMLS.w0,dof_sets_mod,T_LIR] = LIR(K_AMLS.w0,M_AMLS.w0,...
        dof_sets_mod,options);
    
    % update transformation matrix
    T_AMLS.w0 = T_AMLS.w0*T_LIR;
end


%% form residual correction transformation
t_start = tic;
T_AMLS.w2 = Psi_hat*(Frs*(Psi_hat'*(M_free*T_AMLS.w0)));

% print timing info to screen
if options.verbose
    fprintf([tbi,'T_AMLS+ assembly time: %6.2f \n'],toc(t_start))
end

% calculate residual correction mass and stiffness components
t_start = tic;

% residual enhanced stiffness terms
K_AMLS.w4 = T_AMLS.w2'*K_free*T_AMLS.w2;
K_AMLS.w2 = 0;

% residual enhanced mass terms
M_AMLS.w4 = T_AMLS.w2'*M_free*T_AMLS.w2;
M_AMLS.w2 = K_AMLS.w4;

% print timing info to screen
if options.verbose
    fprintf([tbi,'Time to apply transformations to K and M: %6.2f \n'],toc(t_start))
end

% Make sure that matrices that should be symmetric are symmetric
% ======================================================================= %
K_AMLS.w0 = (1/2)*(K_AMLS.w0 + K_AMLS.w0);
K_AMLS.w4 = (1/2)*(K_AMLS.w4 + K_AMLS.w4);

M_AMLS.w0 = (1/2)*(M_AMLS.w0 + M_AMLS.w0);
M_AMLS.w4 = (1/2)*(M_AMLS.w4 + M_AMLS.w4);


% place AMLS transformation in output cell
varargout{1} = T_AMLS;

if options.verbose
    fprintf([tbi,'Residual Enhanced AMLS calculation time: %5.2f s\n'],toc(t_start_AMLSpl))
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
    fprintf([tb,'\n'])
end