% Erreur entre les deux versions en fonction du conditionnement.
% 100 matrices 100x100 tirees aleatoirement, conservees seulement si cond > 1e6.

n = 100;
c = zeros(n, 1);
s = zeros(n, 1);

for i = 1:n
    a = rand(100);
    while cond(a) < 1e6
        a = rand(100);
    end
    x1 = Schmidt_orthogonalization_Original(a);
    x2 = Schmidt_orthogonalization_Modified(a);
    s(i) = Erreur_moyenne(x1, x2);
    c(i) = cond(a);
end

plot(c, s, '--rs', 'LineWidth', 2, ...
     'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'g', ...
     'MarkerSize', 10);
xlabel('condition');
ylabel('erreur');
