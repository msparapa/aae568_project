function [sol] = bvpmc(odefun, quadsfun, bcfun, solinit, options)

sol = solinit;

nOdes = length(sol.y(:,1));

if isempty(quadsfun)
    nQuads = 0;
else
    nQuads = length(sol.quads(:,1));
end

nParams = length(sol.parameters);


options.nOdes = nOdes;
options.nQuads = nQuads;
options.nParams = nParams;
options.odefun = odefun;
options.bcfun = bcfun;
options.consts = sol.consts;
N = options.Nodes;

constraints = @(X)(collocation_constraints(X,options));

x0 = [];
for ii = 1:nOdes
    x0 = [x0, sol.y(ii,:)];
end
x0 = [x0, sol.parameters];

[X,fval,exitflag,output,lambda,grad,hessian] = fmincon(@collocation_cost, x0, [], [], [], [], [], [], constraints, options);

for ii = 1:options.nOdes
    y(ii,:) = X(((ii-1)*N+1):((ii)*N));
end

params = X(options.nOdes*N+1:end);

sol.y = y;
sol.parameters = params;

end

function [c, ceq] = collocation_constraints(X, options)
c = [];
ceq = [];
N = options.Nodes;

for ii = 1:options.nOdes
    y(ii,:) = X(((ii-1)*N+1):((ii)*N));
end

params = X(options.nOdes*N+1:end);

dX = options.odefun(0, y, params, options.consts);

dp0 = dX(:,1:end-1);
dp1 = dX(:,2:end);
p0 = y(:,1:end-1);
p1 = y(:,2:end);

midpoint_predicted = 1/2*(p0+p1) + 1/(N-1)*(dp0-dp1)/8;
midpoint_derivative_predicted = -3/2*(N-1)*(p0-p1) - 1/4*(dp0+dp1);
midpoint_derivative_actual = options.odefun(0, midpoint_predicted, params, options.consts);

collo_constraint = midpoint_derivative_predicted - midpoint_derivative_actual;
ceq = [ceq; collo_constraint(:)];

bcs = options.bcfun(0, y(:,1), 1, y(:,end), 0, 0, params, options.consts);

ceq = [ceq;bcs'];

end

function [J] = collocation_cost(X)
J=0;
end