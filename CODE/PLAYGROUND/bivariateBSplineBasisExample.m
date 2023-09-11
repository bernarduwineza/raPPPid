% Simple example for bivariate B-splines 

uKnot = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
vKnot = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
p     = 3;
q     = 3;
weights = [1 1 1 1 1 1 1];

[X,Y] = meshgrid (linspace(0,1,60));
R     = zeros(size(X,1),size(X,1));
Ki = 3;
for i=1:size(X,1)
    for j=1:size(X,1)
        Xi  = X(1,i);
        Eta = Y(j,1);
        % Ki = FindSpans(3, p, Xi, vKnot);
        Ni = Der1BasisFuns (Ki, Xi, p, vKnot);
        % Ki = FindSpans(3, q, Eta, uKnot);
        Mi = Der1BasisFuns (Ki, Eta, q, uKnot);
        R(i,j) = Ni'*Mi;
        fprintf('i = %d, j= %d \n', i,j);
    end
end

noPts   = 120;
xi      = linspace(0,1,noPts);
N       = zeros(noPts,length(Ni));

for i=1:noPts
    % Ki = FindSpans(3, p, xi(i), vKnot);
    N(i,:) = Der1BasisFuns (Ki, xi(i), p, uKnot);
end

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);

figure
hold on
plot(xi,N(:,1:4),'Linewidth',2.5)

%exportfig(gcf,'cubic-N3',opts)

figure
surf (X,Y,R)
title('Basis function associated to a local knot vector')
hold off