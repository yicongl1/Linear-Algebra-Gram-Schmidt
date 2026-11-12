% Schmidt_orthogonalization_Original.m
function b = Schmidt_orthogonalization_Original(a)
% Orthogonalisation de Schmidt, version classique (CGS).
% a : matrice dont les colonnes sont orthonormalisees
% b : matrice orthonormee de meme taille que a
[m, n] = size(a);
b = zeros(m, n);
b(:,1) = a(:,1);
for i = 2:n
    b(:,i) = a(:,i);
    for j = 1:i-1
        b(:,i) = b(:,i) - dot(a(:,i), b(:,j)) / dot(b(:,j), b(:,j)) * b(:,j);
    end
end
for k = 1:n
    b(:,k) = b(:,k) / norm(b(:,k));
end
end
