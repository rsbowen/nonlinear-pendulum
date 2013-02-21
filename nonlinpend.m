function nonlinpend()
%{
thetainit = @(xs) 0.7*ones(size(xs));
[xs theta] = nonlinonce(thetainit,100,2*pi);
figure(1);
plot(xs,theta)
xlabel('Starting with theta=0.7')

thetainit = @(xs) 0.7*ones(size(xs)) + 3*sin(xs/2);
[xs theta] = nonlinonce(thetainit,100,2*pi);
figure(2);
plot(xs,theta);
xlabel('Starting with theta=0.7+3*sin(theta/2)')

thetainit = @(xs) 0.7*ones(size(xs)) -sin(7*xs/2)
[xs theta] = nonlinonce(thetainit,100,2*pi);
figure(3)
plot(xs,theta);
xlabel('Starting with theta=0.7 - sin(7 theta/2)')
%}
thetainit = @(xs) 0.7*ones(size(xs)) + (pi-0.7)*sin(pi*xs/20);
[xs theta] = nonlinonce(thetainit,100,20,100);
figure(4);
plot(xs,theta);

xlabel('Large T')
%}

end


function [xs theta] = nonlinonce(thetainit,m,T,maxiter)
if(nargin<4)
maxiter =20
end
xs = linspace(0,T,m)';
theta = thetainit(xs);
h = xs(2)-xs(1);
for iter=1:maxiter
    theta = theta+newton(theta,h);
end
end


function delta = newton(theta,h)
J = zeros(numel(theta));
J = diag(-2+h^2 * cos(theta)) + diag(ones(numel(theta)-1,1),1) +...
diag(ones(numel(theta)-1,1),-1);
J = J/(h^2);
%boundary condition
%J = [J;1 zeros(1,numel(theta)-1);zeros(1,numel(theta)-1) 1];
RHS = -(1/h^2)*conv([0.7;theta;0.7],[1 -2 1]);
RHS = RHS(3:end-2) - sin(theta);
%boundary condition
%RHS = [RHS;0;0];
delta = J\RHS;
end


