function [K_CB,M_CB,dof_sets_mod,varargout] = CB(K_free,M_free,dof_sets,varargin)
%
% "Craig-Bampton (CB) Model"
% ==========================
% Dimitri Krattiger (4-3-2017)
%
% description
% ===========
% This code produces reduced mass and stiffness matrices using the 
% Hurty/Craig-Bampton method
%
% inputs
% ======
% K_free    = free stiffness matrix
% 
% M_free    = free mass matrix
%
% i_i       = indices of interior DOFs
%
% i_b       = indices of boundary (interface) DOFs
%
% options   = options structure (must specify an integer value for
%             the number of fixed interface modes, options.n_FI. 
%             Other options are optional)
%
% outputs
% =======
% K_CB      = CB stiffness matrix 
% 
% M_CB      = CB mass matrix 
% 
% T_CB      = Transformation from full DOF vector to CB DOF vector
%             (only computed if options.outputT = true)
%
% Mc        = interior-boundary (ib) partition of CB matrices
%             (only computed if options.resPlus = true)
%
% PHI_FI    = fixed-interface modes (eigenvectors)
%             (only computed if options.resPlus = true)
%
% L_FI      = fixed-interface eigenvalues
%             (only computed if options.resPlus = true)


%% Check what inputs options are given and set rest to default
% ======================================================================= %
% set_default values for options
defaults.n_FI           = [];
defaults.w_i            = [];
defaults.verbose        = false;
defaults.verboseTab     = '%%   ';
defaults.verboseTabSum  = '';
defaults.resPlus        = false;
defaults.plots          = false;
defaults.outputT        = nargout>=3;
defaults.quasiStatic    = false;
defaults.wQS            = 0; %either wQS=0, or 3 frequencies (top-ish, center and bottom-ish)

if ~isstruct(varargin{1})
    options.n_FI = varargin{1};
else
    options = varargin{1};
end

options = setstructfields(defaults,options);

% check if property arguments are acceptable
%CheckInputs(options);

%% Display Timing info
% ======================================================================= %

if options.verbose
    t_start_CB = tic;
    tb = options.verboseTabSum;
    options.verboseTabSum = [options.verboseTabSum,options.verboseTab];
    tbi = options.verboseTabSum;    
    fprintf([tb,'\n'])
    fprintf([tb,'Craig Bampton (CB) calculation\n'])
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])    
end

%% extract interior and boundary dof index sets from dof_sets structure
% ======================================================================= %

i_i = dof_sets.i;
i_b = [];
dof_set_fields = fieldnames(dof_sets);
for i = 1:length(dof_set_fields)
    if ~isequal(dof_set_fields{i},'i');
        i_b = [i_b;dof_sets.(dof_set_fields{i})];
    end
end

% take index of all interior DOFs from "dof_sets" structure
n_i = length(i_i);
n_b = length(i_b);

%% Perform CB reduction
% ======================================================================= %

% extract blocks of mass and stiffness matrices
K_ii = K_free(i_i,i_i);
K_ib = K_free(i_i,i_b);
K_bb = K_free(i_b,i_b);

M_ii = M_free(i_i,i_i);
M_ib = M_free(i_i,i_b);
M_bb = M_free(i_b,i_b);

%% compute fixed interface modes
% ======================================================================= %
% % do full eigenvalue calculation?
do_full_eig = false;
if n_i<800
    do_full_eig = true;
end

% use cutoff frequency?
use_w_cut = true;
if isempty(options.w_i)
    use_w_cut = false;
end

% centering frequency (about which to compute eigenvalues for interior);
wQS = options.wQS;
wCenter = (1/2)*(max(wQS)+min(wQS));
t_start = tic;
if do_full_eig
    
    % if substructure is small enough use direct eigenvalue solver
    [PHI_FI,L_FI] = eig(full(K_ii),full(M_ii));

    
    if options.verbose
        fprintf([tbi,'Full eigenvalue calculation time: %5.2f s\n'],toc(t_start))
    end
    
else
    if use_w_cut
        n_FI_temp = 50;
        eig_opts.p = 3*n_FI_temp;
        [PHI_FI,L_FI] = eigs(K_ii,M_ii,n_FI_temp,'sm',eig_opts);

        % compute more fixed-interface modes if necessary
        while max(abs(sqrt(diag(L_FI))))<options.w_i
            
            % if current fixed interface modes don't reach the frequency
            % cutoff, double the number that are computed
            n_FI_temp = n_FI_temp*2;
            if n_FI_temp > n_i
                warning('no truncation: cutoff frequency, "w_i" is too high for this model. Proceeding Anyway...')
                [PHI_FI,L_FI] = eig(full(K_ii),full(M_ii));
                break
            else
                [PHI_FI,L_FI] = eigs(K_ii,M_ii,n_FI_temp,'sm');
            end
        end
        
    else % no cutoff frequency is available
        n_FI = options.n_FI;
        if n_FI==0
            eig_opts.dummy = [];
        else
            eig_opts.p = min(3*n_FI,size(K_ii,1)); % number of Lanczos vectors
        end
        [PHI_FI,L_FI] = eigs(K_ii,M_ii,n_FI,wCenter^2,eig_opts); 
    end
    
    if options.verbose
        fprintf([tbi,'Iterative eigenvalue calculation time: %5.2f s\n'],toc(t_start))
    end
end

% sort eigenvalues and eigenvectors by distance from centering frequency
L_FI = diag(L_FI);
[~,i_FI] = sort(abs(sqrt(L_FI)-wCenter));
L_FI = L_FI(i_FI);
PHI_FI = PHI_FI(:,i_FI);

% determine number of fixed interface modes
if use_w_cut
    % truncate by frequency cutoff
    n_FI = sum(abs(sqrt((L_FI)))<options.w_i);
    if n_FI==0
        warning('Cutoff frequency may be too low. No fixed interface modes were kept in model')
    end
else
    % truncate by number given in options
    n_FI = options.n_FI;
end

% truncate fixed-interface modes and frequencies to n_FI
L_FI = L_FI(1:n_FI);
PHI_FI = PHI_FI(:,1:n_FI);

% do a normal sort of fixed-interface modes (ignoring centering freqency)
[L_FI,i_FI] = sort(L_FI);
PHI_FI = PHI_FI(:,i_FI);

if options.verbose
    fprintf([tbi,'Kept %i fixed-interface modes in CB model.\n'],n_FI)
end

% mass normalize fixed interface modes (not always necessary but may
% improve conditioning)
PHI_FI = PHI_FI*diag(diag(PHI_FI'*M_ii*PHI_FI).^(-0.5));

t_start = tic;



% fill in wQS values that were assigned as nan
Lmin = min(L_FI);   Lmax  = max(L_FI);
n_wQS = length(wQS);

QS_Sliders = linspace(0.05,0.95,n_wQS).^1;
% QS_Sliders = linspace(-0.05,1.05,n_wQS).^1;
LQS2 = Lmin + (Lmax-Lmin)*QS_Sliders;

% overwrite any values of wQS that are nans
LQS = wQS.^2;
LQS(isnan(LQS)) = LQS2(isnan(LQS));
LQS = sort(LQS);

% round to 5 significant figures and then find unique values
sigfigs = 4;
[~,i_unique] = unique(sd_round(LQS,sigfigs));
LQS = LQS(i_unique);
wQS = sqrt(LQS);
n_wQS = length(LQS);

% compute constraint modes (or quasi static constraint modes)
% Psi = -K_ii\Kib; %if wQS=0
Psi = zeros(n_i,n_b*n_wQS);
for i = 1:n_wQS
    t_startQS = tic;
    D_ii = K_ii-LQS(i)*M_ii;
    D_ib = K_ib-LQS(i)*M_ib;
    
    % possible bug in backslash operator for 
    % symmetric sparse matrices? but mldivide
    % works fine
    Psi(:,(i-1)*n_b + (1:n_b)) = full(-mldivide(D_ii,D_ib)); 
    if options.verbose
        if all(LQS == 0)
            fprintf([tbi,'Constraint Mode Calc. Time: %5.2f s\n'],i,toc(t_startQS))
        else
            fprintf([tbi,'Quasi-Static Constraint Mode Set %i Calc. Time: %5.2f s, freq: %4.2f rad/s\n'],i,toc(t_startQS),sqrt(LQS(i)))
        end
    end
end

if options.verbose
    fprintf([tbi,'Constraint mode calculation time: %5.2f s\n'],toc(t_start))
end

%% Form CB mass and stiffness matrices
% ======================================================================= %

% start timer for M_CB
t_start = tic;

% intermediate steps in forming M_CB
MiiPsi = M_ii*Psi;
MbiPsi = repmat(M_ib,[1,n_wQS])'*Psi;
Mc = (MiiPsi + repmat(M_ib,[1,n_wQS]));
M_CB_ib = PHI_FI'*Mc;
M_CB_bb = Psi'*MiiPsi + MbiPsi + MbiPsi'+repmat(M_bb,[n_wQS,n_wQS]);

% Assemble M_CB
M_CB = [eye(n_FI),  M_CB_ib;...
        M_CB_ib',   M_CB_bb];
 
if options.verbose
    fprintf([tbi,'M_CB assembly time: %5.2f s\n'],toc(t_start))
end

% start timer for K_CB
t_start = tic; 
if isequal(wQS,0)
    
    % Intermediate steps for K_CB
    K_CB_bb = K_ib'*Psi+K_bb;

    % Assemble K_CB
    K_CB = [diag(L_FI(1:n_FI)),zeros(n_FI,n_b);...
            zeros(n_b,n_FI),K_CB_bb];
else
    
    % Intermediate steps for K_CB
    KiiPsi = K_ii*Psi;
    KbiPsi = repmat(K_ib,[1,n_wQS])'*Psi;
    Kc = (KiiPsi + repmat(K_ib,[1,n_wQS]));
    K_CB_ib = PHI_FI'*Kc;
    K_CB_bb = Psi'*KiiPsi + KbiPsi + KbiPsi'+repmat(K_bb,[n_wQS,n_wQS]);

    % Assemble K_CB
    K_CB = [diag(L_FI(1:n_FI)), K_CB_ib;...
            K_CB_ib',           K_CB_bb];
end
    
if options.verbose
    fprintf([tbi,'K_CB assembly time: %5.2f s\n'],toc(t_start))
end

% symmetrize Craig-Bampton mass and stiffness matrices
K_CB = (1/2)*(K_CB + K_CB');
M_CB = (1/2)*(M_CB + M_CB');

%% Form updated DOF set structure
% ======================================================================= %

% internal dofs replaced with fixed interface modal DOFs


% generate dof_set structure with empty fields
for i = 1:length(dof_set_fields)    
    dof_sets_mod.(dof_set_fields{i}) = [];
end
dof_sets_mod.i = (1:n_FI)';

% boundary DOFs are grouped by sorted by  quasi-static frequency and then
% by boundary set
count = n_FI;

% need to assign the entire dof set structure on
for j = 1:n_wQS
    for i = 1:length(dof_set_fields)
        if ~isequal(dof_set_fields{i},'i');
            n_dof_set = length(dof_sets.(dof_set_fields{i}));
            dof_sets_mod.(dof_set_fields{i}) = ...
                [dof_sets_mod.(dof_set_fields{i}); count + (1:n_dof_set)'];
            count = count + n_dof_set;
        end
    end
end

%% Form BMS transformation (CB transformation)
% ======================================================================= %
if options.outputT
    
    % Assemble CB transformation by concatenating appropriate blocks
    T_CB = [PHI_FI,Psi;zeros(n_b,n_FI),repmat(eye(n_b),[1,n_wQS])];

    % re-sort transformation rows to original DOF sort
    % (must do this so that rows of transformation matrix correspond to rows of
    %  original mass and stiffness matrices)
    [~,i_shuffle1] = sort([i_i;i_b]);
    T_CB = T_CB(i_shuffle1,:);
    
    % store in variable size output array
    varargout{1} = T_CB;
else
    varargout{1} = [];
end

% optional outputs that are necessary to perform residual enhancement
if options.resPlus
    
    % interior-boundary (ib) partition of CB mass
    varargout{2} = Mc;
    
    % fixed-interface modes (eigenvectors)
    varargout{3} = PHI_FI;
    
    % fixed-interface eigenvalues
    varargout{4} = L_FI;
    
    % interior-boundary (ib) partition of CB stiffness
    % (only non-zero if quasi-static constraint modes are used)
    if ~isequal(wQS,0)    
        varargout{5} = Kc;
    end
else
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
    varargout{5} = [];    
end

if nargout>=9
    varargout{6} = wQS;
end

%% Display Timing info
% ======================================================================= %

if options.verbose
    fprintf([tbi,'Craig-Bampton calculation time: %5.2f s\n'],toc(t_start_CB))
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
    fprintf([tb,'\n'])
end