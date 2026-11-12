% Comparaison des deux versions sur les matrices du rapport.
% Resultats attendus :
%   cond(A) = 53.9413,  e(A) ~ 0
%   cond(B) = 1.7321e+10,  e(B) = 0.2679

e = 1e-10;

A = [1 2 3; 2 -1 2; 1 1 2];
A1 = Schmidt_orthogonalization_Original(A);
A2 = Schmidt_orthogonalization_Modified(A);
eA = Erreur_moyenne(A1, A2);

B = [1 1 1; e 0 0; 0 e 0; 0 0 e];
B1 = Schmidt_orthogonalization_Original(B);
B2 = Schmidt_orthogonalization_Modified(B);
eB = Erreur_moyenne(B1, B2);

fprintf('A1 (version originale) :\n');
disp(A1);
fprintf('A2 (version modifiee) :\n');
disp(A2);
fprintf('cond(A) = %.4f\n', cond(A));
fprintf('e(A)    = %.4e\n\n', eA);

fprintf('B1 (version originale) :\n');
disp(B1);
fprintf('B2 (version modifiee) :\n');
disp(B2);
fprintf('cond(B) = %.4e\n', cond(B));
fprintf('e(B)    = %.4f\n', eB);
