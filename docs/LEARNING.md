# 算法背景与代码笔记 · Contexte algorithmique et notes de code

[Français](#notes-de-code) · [中文](#中文笔记)

本文档记录项目的数学背景、相关理论，以及两种实现的代码要点。完整推导与实验数据见
[`../report/rapport-TP-2022.pdf`](../report/rapport-TP-2022.pdf)。

Ce document rassemble le contexte mathématique, les notions associées et les points clés des deux
implémentations. La théorie complète et les résultats expérimentaux figurent dans le
[rapport](../report/rapport-TP-2022.pdf).

---

## 中文笔记

### 1. 问题背景

给定欧氏空间中的线性无关族 `(v₁,…,vₙ)`，正交化的目标是构造标准正交族 `(e₁,…,eₙ)`，使得

$$\operatorname{Vect}(v_1,\dots,v_k) = \operatorname{Vect}(e_1,\dots,e_k), \quad k = 1,\dots,n.$$

每一步的核心运算是向已有向量的投影：

$$\operatorname{proj}_{u}(v) = \frac{\langle v, u\rangle}{\langle u, u\rangle}\,u.$$

两个版本的差别只在投影的组织方式：

- **经典版（CGS）**：先完整计算 `u_k = v_k − Σ_{j<k} ⟨v_k, u_j⟩ u_j`，再归一化；
- **修正版（MGS）**：按 `j = 1, 2, …` 的顺序，每次投影后立刻用更新过的向量继续下一次投影。

两者在精确算术下等价；在浮点运算下，CGS 的误差会沿着迭代链传递并放大，MGS 则把误差限制
在单次投影内。

### 2. 相关理论

- **秩定理**：`dim ker A + rang A = 列数`，用于判断一族向量的线性相关性；
- **投影与正交补**：`ker A ⊕ (ker A)^⊥ = ℝⁿ`，正交化本质上是逐步构造正交补中的分量；
- **条件数**：`cond(A) = ‖A‖·‖A⁻¹‖`，衡量线性问题对扰动的敏感程度。条件数越大，列向量越
  接近线性相关，投影中的相消越严重；
- **Frobenius 范数**：`‖E‖_F = √(Σ e_ij²)`，用于度量两种算法输出矩阵之间的差异；
- **实对称矩阵**：可正交对角化 `A = PDPᵀ`，其特征子空间两两正交，是 Gram–Schmidt 的典型
  应用场景。

### 3. C 实现要点

`main.c` 由两个辅助函数和一个主循环组成。

```c
double produitscalaire(double a[], double b[], int l)  // 内积 ⟨a,b⟩
{
    double s = 0;
    for (int i = 1; i <= l; i++) s += a[i] * b[i];
    return s;
}

void norme1(double a[], int l)                         // 单位化 a/‖a‖
{
    double s = 0;
    for (int i = 1; i <= l; i++) s += a[i] * a[i];
    s = sqrt(s);
    for (int i = 1; i <= l; i++) a[i] /= s;
}
```

主循环对应 CGS 的公式：

```c
for (int i = 1; i <= y; i++) {
    for (int j = 1; j <= x; j++) b[i][j] = a[i][j];       // u_k ← v_k
    for (int j = 1; j < i; j++) {
        double ps = produitscalaire(a[i], b[j], x);       // ⟨v_k, u_j⟩
        for (int k = 1; k <= x; k++)
            b[i][k] -= ps * b[j][k];                      // 减去投影分量
    }
    norme1(b[i], x);
}
```

实现上的几个选择：

- **转置读入**：输入时写入 `a[j][i]` 而非 `a[i][j]`，使原矩阵的每一列在数组中成为连续的一行，
  后续按 `a[i]` 取向量即为按列处理；
- **1-based 下标**：数组声明为 `100×100`，仅使用 `1..n`，让代码下标与数学记号一致；
- **边界**：固定上限 100，仅适用于方阵；`double a[100][100] = {}` 依赖 C23 / GNU 扩展，
  可移植写法为 `= {{0}}`。

编译与测试：

```bash
gcc -O2 -Wall -o gs main.c -lm
printf "3 3\n1 2 3\n2 -1 2\n1 1 2\n" | ./gs
```

输出与 MATLAB 的 `Q` 完全一致，说明 C 版本实现了正确的正交化。

### 4. MATLAB 实现要点

完整文件位于 [`../matlab/`](../matlab/)：
CGS 与 MGS 各一个函数文件、误差度量 `Erreur_moyenne.m`，以及 `demo_comparison.m`、
`condition_experiment.m` 两个脚本。

误差度量函数：

```matlab
function s = Erreur_moyenne(a, b)
    d = a - b;
    s = 0;
    [m, n] = size(a);
    for i = 1:m
        for j = 1:n
            s = s + d(i, j)^2;
        end
    end
end
```

该函数实际返回 `‖a − b‖_F²`，缺少最后一步开方。由于平方根单调，误差随条件数的变化趋势与
结论不受影响；报告中将其简称为 Frobenius 范数。

条件数实验：

```matlab
c = zeros(100, 1);  s = zeros(100, 1);
for i = 1:100
    a = rand(100);
    while cond(a) < 1e6            % 只保留病态矩阵
        a = rand(100);
    end
    x1 = Schmidt_orthogonalization_Original(a);
    x2 = Schmidt_orthogonalization_Modified(a);
    s(i) = Erreur_moyenne(x1, x2);
    c(i) = cond(a);
end
plot(c, s, '--rs', 'LineWidth', 2, ...
     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10)
xlabel('condition'); ylabel('erreur');
```

实验结果显示，`cond` 增大时两种版本的差异随之增大，实验图见
[`../report/figures/`](../report/figures/)。

---

## Notes de code

### 1. Contexte du problème

Étant donnée une famille libre `(v₁,…,vₙ)` d'un espace euclidien, l'orthogonalisation construit une
famille orthonormée `(e₁,…,eₙ)` telle que

$$\operatorname{Vect}(v_1,\dots,v_k) = \operatorname{Vect}(e_1,\dots,e_k), \quad k = 1,\dots,n.$$

L'opération de base est la projection sur un vecteur déjà construit :

$$\operatorname{proj}_{u}(v) = \frac{\langle v, u\rangle}{\langle u, u\rangle}\,u.$$

Les deux versions ne diffèrent que par l'organisation de ces projections :

- **version classique (CGS)** : calcul complet de `u_k = v_k − Σ_{j<k} ⟨v_k, u_j⟩ u_j`, puis
  normalisation ;
- **version modifiée (MGS)** : pour `j = 1, 2, …`, chaque projection est suivie immédiatement de la
  suivante, sur le vecteur déjà mis à jour.

Les deux sont équivalentes en arithmétique exacte ; en arithmétique flottante, l'erreur de CGS se
propage le long de la chaîne d'itérations, tandis que celle de MGS reste locale à chaque projection.

### 2. Notions associées

- **Théorème du rang** : `dim ker A + rang A = nombre de colonnes`, pour tester la dépendance
  linéaire d'une famille ;
- **Projection et orthogonal** : `ker A ⊕ (ker A)^⊥ = ℝⁿ` ;
- **Conditionnement** : `cond(A) = ‖A‖·‖A⁻¹‖`, mesure de la sensibilité aux perturbations. Plus il
  est grand, plus les colonnes sont proches de la dépendance linéaire et plus les cancellations sont
  importantes ;
- **Norme de Frobenius** : `‖E‖_F = √(Σ e_ij²)`, pour mesurer l'écart entre deux matrices de sortie ;
- **Matrices symétriques réelles** : diagonalisables dans une base orthonormée (`A = PDPᵀ`), cas
  d'application typique de Gram-Schmidt.

### 3. Points clés de l'implémentation C

`main.c` se compose de deux fonctions auxiliaires et d'une boucle
principale.

```c
double produitscalaire(double a[], double b[], int l)  /* ⟨a,b⟩ */
{
    double s = 0;
    for (int i = 1; i <= l; i++) s += a[i] * b[i];
    return s;
}

void norme1(double a[], int l)                         /* a / ‖a‖ */
{
    double s = 0;
    for (int i = 1; i <= l; i++) s += a[i] * a[i];
    s = sqrt(s);
    for (int i = 1; i <= l; i++) a[i] /= s;
}
```

La boucle principale traduit la formule de CGS :

```c
for (int i = 1; i <= y; i++) {
    for (int j = 1; j <= x; j++) b[i][j] = a[i][j];       /* u_k ← v_k */
    for (int j = 1; j < i; j++) {
        double ps = produitscalaire(a[i], b[j], x);       /* ⟨v_k, u_j⟩ */
        for (int k = 1; k <= x; k++)
            b[i][k] -= ps * b[j][k];                      /* retrait de la projection */
    }
    norme1(b[i], x);
}
```

Choix d'implémentation :

- **lecture transposée** : l'entrée est stockée dans `a[j][i]`, de sorte que chaque colonne de la
  matrice d'origine devienne une ligne contiguë du tableau ;
- **indices à partir de 1** : tableaux de taille `100×100`, seule la partie `1..n` est utilisée, ce
  qui aligne le code sur les notations mathématiques ;
- **limites** : taille fixe, matrices carrées uniquement ; `double a[100][100] = {}` relève de C23
  ou d'une extension GNU, la forme portable étant `= {{0}}`.

```bash
gcc -O2 -Wall -o gs main.c -lm
printf "3 3\n1 2 3\n2 -1 2\n1 1 2\n" | ./gs
```

La sortie coïncide avec la matrice `Q` obtenue sous MATLAB.

### 4. Points clés de l'implémentation MATLAB

Les fichiers complets se trouvent dans [`../matlab/`](../matlab/) :
une fonction par version, la mesure d'erreur `Erreur_moyenne.m`, puis les scripts
`demo_comparison.m` et `condition_experiment.m`.

Mesure d'erreur :

```matlab
function s = Erreur_moyenne(a, b)
    d = a - b;
    s = 0;
    [m, n] = size(a);
    for i = 1:m
        for j = 1:n
            s = s + d(i, j)^2;
        end
    end
end
```

La fonction renvoie en réalité `‖a − b‖_F²`, sans la racine carrée finale. La racine étant monotone,
l'évolution de l'erreur avec le conditionnement et les conclusions restent inchangées ; le rapport
la désigne par simplification comme norme de Frobenius.

```matlab
for i = 1:100
    a = rand(100);
    while cond(a) < 1e6            % conservation des seules matrices mal conditionnées
        a = rand(100);
    end
    x1 = Schmidt_orthogonalization_Original(a);
    x2 = Schmidt_orthogonalization_Modified(a);
    s(i) = Erreur_moyenne(x1, x2);
    c(i) = cond(a);
end
plot(c, s, '--rs', 'LineWidth', 2, ...
     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10)
xlabel('condition'); ylabel('erreur');
```

L'écart entre les deux versions croît avec `cond` ; les figures correspondantes se trouvent dans
[`../report/figures/`](../report/figures/).
