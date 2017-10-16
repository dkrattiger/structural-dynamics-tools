function [K_CB,M_CB,dof_sets_mod,T_CB] = CB_plus_LIR(K, M, dof_sets, options)

% "Craig-Bampton (CB) with residual mode enhancement (plus) and 
% Local Interface Reduction (LIR)"
% ==========================
% Dimitri Krattiger (4-3-2017)
%
% description
% ===========
% This code produces reduced mass and stiffness matrices using the enhanced
% CB method discussed in the reference:
%
%       Kim, J.G. and Lee, P.S. An Enhanced Craig-Bampton Method -
%       International Journal for Numerical Methods in Engineering, 2015
%
% inputs
% ======
% K         = stiffness matrix
% 
% M         = mass matrix
%
% dof_sets  = structure containing indices of node sets (obtain with 
%             "find_node_sets.m" and "node2dof.m" functions)
% 
% options   = options structure. 
%             Must specify the number of fixed interface modes or the 
%             interior cutoff frequency  (options.n_FI -or- options.w_i)
%             Must specify the number of boundary modes or the 
%             boundary cutoff frequency  (options.n_CC -or- options.w_b)
%            
%
% outputs
% =======
% K_CB      = CB stiffness matrix with interface reduction and residual 
%             enhancement
% 
% M_CB      = CB mass matrix with interface reduction and residual 
%             enhancement
%
% dof_sets_mod  = new dof sets structure for CB matrices
% 
% T_CB      = Transformation from full DOF vector to CB DOF vector

%% Check what inputs are given and set rest to default
% ======================================================================= %

% set_default values for options
defaults.BoundaryMethod = 'none';
defaults.n_FI           = [];
defaults.n_CC           = [];
defaults.w_i            = [];
defaults.w_b            = [];
defaults.verbose        = false;
defaults.verboseTab     = '%%   ';
defaults.verboseTabSum  = '';
defaults.plots          = false;
defaults.outputT        = nargout>=4;
defaults.wQS            = 0;

% overwrite defaults with any options that are input and copy defaults
% where no options are specified
options = setstructfields(defaults,options);

%% Display Timing info
% ======================================================================= %

if options.verbose
    t_start_CBpl = tic;    
    tb = options.verboseTabSum;
    options.verboseTabSum = [options.verboseTabSum,options.verboseTab];
    tbi = options.verboseTabSum;
    fprintf([tb,'\n'])
    fprintf([tb,'Residual Enhanced CB (CB+) Reduction\n'])
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
end

%% perform basic CB transformation
% ======================================================================= %

% extract interior and boundary dof sets
i_i = dof_sets.i;
dof_set_names = fieldnames(dof_sets);
i_b = [];
for i = 2:length(dof_set_names)
    i_b = [i_b;dof_sets.(dof_set_names{i})];
end

% Original number of boundary and interface DOFs
n_b = length(i_b);
n_i = length(i_i);

% CB reduction options
options_CB              = options;
options_CB.outputT      = true;
options_CB.resPlus      = true;

if isequal(options.wQS,0)
    [K_CB.w0,M_CB.w0,dof_sets_mod,T_CB.w0,Mc,PHI_d,L_d] = CB(K,M,dof_sets,options_CB);
else
    [K_CB.w0,M_CB.w0,dof_sets_mod,T_CB.w0,Mc,PHI_d,L_d,Kc] = CB(K,M,dof_sets,options_CB);    
end

% number of fixed interface modes
n_FI = length(dof_sets_mod.i);

%% Perform LIR reduction
% ======================================================================= %
if ~strcmpi(options.BoundaryMethod,'none')
    
    % options for LIR reduction
    options_LIR             =  options;
    options_LIR.outputT     =  true;
    
    % perform LIR reduction
    [K_CB.w0,M_CB.w0,dof_sets_mod,T_LIR] = LIR(K_CB.w0,M_CB.w0,...
        dof_sets_mod,options_LIR);
    
    % update transformation matrix
    T_CB.w0 = T_CB.w0*T_LIR;    
    Mc = Mc*T_LIR((n_FI+1):end,(n_FI+1):end);
    if ~isequal(options.wQS,0)
        Kc = Kc*T_LIR((n_FI+1):end,(n_FI+1):end);
    end
end

% get new boundary size (whether or not it has changed)
n_b2 = size(K_CB.w0,2)-n_FI;

%% Enhance CB-LIR reduced matrices
% ======================================================================= %

% Could perform the CB enhancement in the CB code, however these
% computations are much faster after LIR has been performed.
t_start = tic;
    
% boundary dof index in CB reduced model
i_b2 = (n_FI+1):size(M_CB.w0,1);
i_FI = 1:n_FI;
M_CB_ib = M_CB.w0(i_FI,i_b2);
K_CB_ib = K_CB.w0(i_FI,i_b2);

% form residual-enrichment transformation (and re-shuffle)
K_ii = K(i_i,i_i);
M_ii = M(i_i,i_i);

% intermediate matrices for performing CB-LIR enhanced matrices
FrsMc = K_ii\Mc - PHI_d*diag(L_d.^-1)*full(M_CB_ib);
McFrsMc = (Mc')*FrsMc;
MiiFrsMc = M_ii*FrsMc;
McFrsMiiFrsMc = FrsMc'*MiiFrsMc;

if ~isequal(options.wQS,0)
    FrsKc = K_ii\Kc - PHI_d*diag(L_d.^-1)*full(K_CB_ib);
end

% form enhanced CB transformation
% if nargout>=4  
if true
    % re-sorting index to properly sort the rows of the transformation
    % matrices
    [~,i_shuffle] = sort([i_i;i_b]);tic
    
    % term to add into T_CB.w0
    if ~isequal(options.wQS,0)
        T_CB_w0_add = [zeros(n_i,n_FI),-FrsKc;zeros(n_b,n_FI + n_b2)];  
        T_CB_w0_add = T_CB_w0_add(i_shuffle,:);
        T_CB.w0 = T_CB.w0+T_CB_w0_add;
    end
    
    % Matrix has zero blocks, but they are small (relatively) so "sparse"
    % is not advantageous
    % T_CB.w2 = sparse([zeros(n_i,n_FI),FrsMc;zeros(n_b,n_FI + n_b2)]);
    T_CB.w2 = [zeros(n_i,n_FI),FrsMc;zeros(n_b,n_FI + n_b2)];    
    
    % re-shuffle the rows of T_CB to match the original input DOF sorting
    [~,i_shuffle] = sort([i_i;i_b]);tic
    T_CB.w2 = T_CB.w2(i_shuffle,:);
end


if isequal(options.wQS,0)
    % partial computation of mass matrix terms
    TcbMTr = [zeros(n_FI,n_FI+n_b2);zeros(n_b2,n_FI),McFrsMc];
    TrMTr = McFrsMiiFrsMc;

    % Residual enhanced mass assembly
    M_CB.w4 = [zeros(n_FI,n_FI+n_b2);zeros(n_b2,n_FI),TrMTr];
    M_CB.w2 = TcbMTr;

    % assemble residual enhanced stiffness terms
    TrKTr = McFrsMc;
    K_CB.w2 = 0;
    K_CB.w4 = [zeros(n_FI,n_FI+n_b2);zeros(n_b2,n_FI),TrKTr];

    % symmetrize Craig-Bampton mass and stiffness matrices
    K_CB.w4 = (1/2)*(K_CB.w4+K_CB.w4');
    M_CB.w4 = (1/2)*(M_CB.w4+M_CB.w4');
else
    
    % SLOW SLOW SLOW SLOW SLOW
    K_CB.w0 = T_CB.w0'*K*T_CB.w0;
    M_CB.w0 = T_CB.w0'*M*T_CB.w0;
    
    K_CB.w2 = T_CB.w0'*K*T_CB.w2;
    M_CB.w2 = T_CB.w0'*M*T_CB.w2;
    
    K_CB.w4 = T_CB.w2'*K*T_CB.w2;
    M_CB.w4 = T_CB.w2'*M*T_CB.w2;
end

% output timing info
if options.verbose
    fprintf([tbi,'M_CB+, K_CB+ calculation time: %5.2f s\n'],toc(t_start))
end


% The calculations in this section are simplified compared to the theory
% presented in "An enhanced Craig–Bampton method" by J.G. Kim and P.S. Lee.
% The matrix simplifications are shown in "A simplified error estimator for 
% the CB method and its application to error control" by S.H. Boo, 
% J.G. Kim, and P.S. Lee.

%% if necessary output transformation matrix
% ======================================================================= %
if ~defaults.outputT
    T_CB = [];
end

if options.verbose
    fprintf([tbi,'Residual Enhanced CB calculation time: %5.2f s\n'],toc(t_start_CBpl))
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
    fprintf([tb,'\n'])
end