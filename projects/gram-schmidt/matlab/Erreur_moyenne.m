function s = Erreur_moyenne(a, b)
% Ecart entre deux matrices, au sens du carre de la norme de Frobenius.
% La racine carree finale est omise, comme dans le rapport : la fonction
% renvoie sum((a-b).^2), ce qui ne change pas l'ordre des resultats.
d = a - b;
s = 0;
[m, n] = size(a);
for i = 1:m
    for j = 1:n
        s = s + d(i,j)^2;
    end
end
end
