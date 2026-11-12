# Gram-Schmidt : classique vs. modifié

**施密特正交化：经典版与修正版**

[Français](#français) · [中文](#中文)

---

## Français

Projet d'algèbre linéaire (2022), réalisé par **LI Yicong** et **CHEN Nuo**.

### Objet

Comparaison de deux algorithmes d'orthonormalisation équivalents en théorie :

- **classique (CGS)** :

  $$u_k = v_k - \sum_{j<k}\langle v_k, u_j\rangle\,u_j$$

- **modifié (MGS)** : le vecteur est mis à jour après chaque projection,

  $$v_k \leftarrow v_k - \langle v_k, u_j\rangle\,u_j, \quad j = 1,\dots,k-1$$

  puis normalisé.

En arithmétique exacte, les deux versions donnent la même base orthonormée `Q`. En arithmétique
flottante, CGS accumule les erreurs d'arrondi d'une étape sur l'autre, tandis que MGS les confine à
chaque projection.

### Résultats

| Matrice | `cond(A)` | `‖Q_CGS − Q_MGS‖_F` |
|---|---|---|
| `A` (3×3, ci-dessous) | 53,9413 | 5,9242×10⁻³¹ |
| `B` (4×3, `ε = 10⁻¹⁰`) | 1,7321×10¹⁰ | 0,2679 |

$$A = \begin{pmatrix}1&2&3\\2&-1&2\\1&1&2\end{pmatrix},\qquad
B = \begin{pmatrix}1&1&1\\ \varepsilon&0&0\\ 0&\varepsilon&0\\ 0&0&\varepsilon\end{pmatrix}$$

Sur 100 matrices aléatoires `100×100` de `cond > 10⁶`, l'écart entre les deux versions croît avec le
conditionnement (`report/figures/fig-erreur-condition-*.png`). La version modifiée est donc
numériquement plus stable.

### Contenu

```
gram-schmidt/
├── src/main.c              # implémentation C (orthonormalisation par colonnes)
├── matlab/                 # implémentations MATLAB et scripts d'expérience
├── matrix.xcodeproj/       # projet Xcode (groupe src)
├── report/
│   ├── rapport-TP-2022.pdf # rapport complet (français)
│   └── figures/            # figures du rapport, sorties MATLAB
└── README.md
```

### MATLAB

```
matlab/
├── Schmidt_orthogonalization_Original.m   # version classique (CGS)
├── Schmidt_orthogonalization_Modified.m   # version modifiée (MGS)
├── Erreur_moyenne.m                       # écart ‖a-b‖_F² entre deux matrices
├── demo_comparison.m                      # reproduit les résultats A et B ci-dessus
└── condition_experiment.m                 # courbes conditionnement — erreur
```

`demo_comparison.m` reproduit `Q_CGS`, `Q_MGS`, `e(A)`, `e(B)` et les conditionnements du
tableau ; `condition_experiment.m` génère les courbes de `report/figures/`.

### Compilation et exécution

```bash
gcc -O2 -Wall -o gs src/main.c -lm
printf "3 3\n1 2 3\n2 -1 2\n1 1 2\n" | ./gs
```

Entrée : nombre de colonnes, nombre de lignes, puis les coefficients ligne par ligne.
L'implémentation actuelle suppose une matrice carrée (`n ≤ 100`) et renvoie la base orthonormée
colonne par colonne :

```
  0.408248  0.759072  0.507093
  0.816497 -0.552052  0.169031
  0.408248  0.345033 -0.845154
```

### Rapport

Le rapport contient la théorie, les deux implémentations MATLAB, les mesures d'erreur et
l'interprétation géométrique (vecteurs presque parallèles, résidus, conditionnement).

---

## 中文

2022 年线性代数课程项目，由**李毅聪**与**陈诺**合作完成。

### 目标

比较两种理论上等价的施密特正交化算法：

- **经典版（CGS）**：一次减去所有已构造向量上的投影；
- **修正版（MGS）**：每做一次投影就立即更新当前向量，再进入下一次投影。

精确算术下两者给出同一个标准正交基 `Q`；浮点运算下，经典版会跨步骤累积舍入误差，修正版
则把误差限制在每次投影内部。

### 结果

| 矩阵 | `cond(A)` | `‖Q_CGS − Q_MGS‖_F` |
|---|---|---|
| `A`（3×3） | 53.9413 | 5.9242×10⁻³¹ |
| `B`（4×3，`ε = 10⁻¹⁰`） | 1.7321×10¹⁰ | 0.2679 |

对 100 个随机生成的 `100×100`、`cond > 10⁶` 的矩阵绘制"条件数—误差"图，误差随条件数
单调上升，说明 MGS 的数值稳定性更好。

### 实现

```
gram-schmidt/
├── src/main.c              # C 实现（按列正交化）
├── matlab/                 # MATLAB 实现与实验脚本
├── matrix.xcodeproj/       # Xcode 工程
├── report/
│   ├── rapport-TP-2022.pdf # 报告全文（法语）
│   └── figures/            # 报告插图与 MATLAB 输出
└── README.md
```

MATLAB 文件：

- `Schmidt_orthogonalization_Original.m`（经典版）、`Schmidt_orthogonalization_Modified.m`（修正版）
- `Erreur_moyenne.m`：两矩阵之差的平方 Frobenius 范数
- `demo_comparison.m`：复现上表的 `Q`、`e(A)`、`e(B)` 与条件数
- `condition_experiment.m`：生成 `report/figures/` 中的条件数—误差曲线

### 编译与运行

```bash
gcc -O2 -Wall -o gs src/main.c -lm
printf "3 3\n1 2 3\n2 -1 2\n1 1 2\n" | ./gs
```

输入为列数、行数，然后逐行输入矩阵；当前实现要求方阵（最大 `100×100`），输出为按列排列的
标准正交基，与报告中的 MATLAB 结果一致。

### 报告

[`report/rapport-TP-2022.pdf`](report/rapport-TP-2022.pdf) 包含理论推导、两版 MATLAB 实现、
误差度量与几何解释。
