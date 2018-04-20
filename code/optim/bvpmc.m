function [sol] = bvpmc(odefun, quadsfun, bcfun, solinit, options)
% Same inputs as bvp4c, but quadsfun = [] for complex (non-reducible)
% systems

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
if options.isdirect
    options.nControls = length(sol.control(:,1));
end
options.consts = sol.consts;
N = options.Nodes;

constraints = @(X)(collocation_constraints(X, options));
cost = @(X)(collocation_cost(X, options));

x0 = [];
for ii = 1:nOdes
    x0 = [x0, sol.y(ii,:)];
end
x0 = [x0, sol.parameters'];

if options.isdirect
    for ii = 1:options.nControls
        x0 = [x0, sol.control(ii,:)];
    end
end

[X,fval,exitflag,output,lambda,grad,hessian] = fmincon(cost, x0, [], [], [], [], [], [], constraints, options);

sol = unwrap_params(X, options);


end

function [X] = wrap_params(sol, options)
% Do I need this?
end

function [sol] = unwrap_params(X, options)
N = options.Nodes;

for ii = 1:options.nOdes
    y(ii,:) = X(((ii-1)*N+1):((ii)*N));
end

params = X(options.nOdes*N+1:(options.nOdes*N)+options.nParams);

if options.isdirect
    b = (options.nOdes*N)+options.nParams;
    for ii = 1:options.nControls
        control(ii,:) = X(b+1+(ii-1)*N:b+(ii)*N);
    end
end
sol.y = y;
sol.parameters = params;
sol.control = control;

end

function [c, ceq] = collocation_constraints(X, options)
c = [];
ceq = [];
N = options.Nodes;

sol = unwrap_params(X, options);
x = linspace(0,1,N);
y = sol.y;
params = sol.parameters;
control = sol.control;

dX = options.odefun(x, y, control, params, options.consts);

dp0 = dX(:,1:end-1);
dp1 = dX(:,2:end);
p0 = y(:,1:end-1);
p1 = y(:,2:end);

x_midpoint = (x(2:end) + x(1:end-1))/2;

midpoint_predicted = 1/2*(p0+p1) + 1/(N-1)*(dp0-dp1)/8;
midpoint_derivative_predicted = -3/2*(N-1)*(p0-p1) - 1/4*(dp0+dp1);
midpoint_derivative_actual = options.odefun(x_midpoint, midpoint_predicted, control(:,1:end-1), params, options.consts);

collo_constraint = midpoint_derivative_predicted - midpoint_derivative_actual;
ceq = [ceq; collo_constraint(:)];

bcs = options.bcfun(x(1), y(:,1), control(:,1), x(end), y(:,end), control(:,end), 0, 0, params, options.consts);

ceq = [ceq;bcs];

end

function [J] = collocation_cost(X, options)
sol = unwrap_params(X, options);
x = linspace(0,1,options.Nodes);
J = options.cost(x, sol.y, sol.parameters, options.consts);
end