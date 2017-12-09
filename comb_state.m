function A = comb_state(n,opt)        


if nargin < 2
    opt = 'num';
end

state = cell(1,n);

if strcmpi(opt,'num')
    for i=1:n
        state{i}=[0 1];
    end
elseif strcmpi(opt,'str')
    for i=1:n
        state{i}='01';
    end
else
    error('Opt for function comb_state can only be num or str.\n');
end

ii = 1:n;
if n==1
   A = state{1}(:);
else
  [A{ii}] = ndgrid(state{ii});
   % concatenate
   A = reshape(cat(n+1,A{:}),[],n) ;
end