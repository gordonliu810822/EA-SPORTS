function options = cvSparseLDASet(opts)

% Set default options.

options.maxIters = 100;
options.tol = 1e-10;
options.epsilon = 0.001;
options.nlam = 20;
options.maxLam = 2;
options.nfold = 5;

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
