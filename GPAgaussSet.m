function options = GPAgaussSet(opts)

%--------------------------------------------------------------------------
% GPAgaussSet creates or alters an options structure for GPAgauss.m.
%--------------------------------------------------------------------------
%   options = GPAgaussSet; (with no input arguments)
%   creates a structure with all fields set to their default values.
%   Each field is an option (also called a parameter).
%
%   GPAgaussSet (with no input or output arguments)
%   displays all options and their default values.
%
%   options = GPAgaussSet(opts); 
%   creates a structure with all fields set to their default values,
%   except valid fields in the structure "opts" replace the defaults.
%
% options.verbose     indicator whether the output is shown in each
%                     iteration, the default is 1.
% options.maxIters    The maximum number of iterations, default is 2000.

% Set default options.

options.verbose = 1;
options.maxIters = 5000;
options.epsZ = 0.005;
options.epsStopLogLik = 1e-10;
options.epsStopPE = 1e-10;
options.initBetaMean = 0.1;
options.initPi = 0.1;
options.initQ1 = 0.75;
options.lbPi1 = 0.001;
options.lbBetaAlpha = 0.001;
options.lbQ1 = 0.001;
options.hom =0;           %indicator whehter the homogeneous variance is used.
options.constraintPi = 0; %indicator whehter the constrain model on pi is used.
options.constraintMu = 0; %indicator whehter the constrain model on mu is used.
options.empiricalNull= 0; %if 1, will use beta(alpha0,1) for null distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End default options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick return if no user opts
  if nargin == 0 || isempty(opts)
    if nargout == 0    % Display options.
      disp('pdco default options:')
      disp( options )
    end
    return
  end

% List of valid field names
  vfields = fieldnames( options );

% Grab valid fields from user's opts
  for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
      options.(field) = opts.(field);
    end
  end
