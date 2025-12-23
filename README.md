# SARAR 模型实证似然方法代码

本目录包含用于 SARAR（空间自回归-自相关）模型的实证似然（Empirical Likelihood, EL）和主成分实证似然（Principal Component Empirical Likelihood, PCEL）方法的 R 代码实现。

## 目录结构

```
PCEL/
├── README.md                    # 本文件
├── code/                        # 主要代码目录
│   ├── RealSimPCEL.R          # PCEL 方法在真实数据上的实现(Sigma估计方式1)
│   ├── RealSimPCEL2.R         # PCEL 方法在真实数据上的实现(Sigma估计方式2)
│   ├── RealSimEL.R            # EL 方法在真实数据上的实现
│   ├── RealSimMAE.R           # 平均绝对误差计算
│   ├── GlambdaChen.R          # Lambda 计算辅助函数
│   ├── PCAEL.R                # PCA 经验似然统计量计算
│   ├── PCA_Data.R             # PCA 数据处理
│   ├── PCA_Data2.R            # PCA 数据处理（第二版）
│   ├── PCA_latex.R            # PCA LaTeX 输出
│   ├── PCA_latex2.R           # PCA LaTeX 输出（第二版）
│   ├── PCA_QQ.R               # PCA Q-Q 图绘制
│   ├── fig4-ab.R              # 图4绘制
│   ├── fig5-abc.R             # 图5绘制
│   ├── PCEL.Rproj             # R 项目文件
│   ├── columbus.dbf           # Columbus 空间数据（DBF格式）
│   ├── columbus.shp           # Columbus 空间数据（SHP格式）
│   └── columbus.shx           # Columbus 空间数据（SHX格式）
├── core/                       # 核心函数目录
│   ├── GlambdaChen.R          # Lambda 计算核心函数
│   └── PCEL.R                 # PCEL 核心函数
├── data/                       # 模拟数据结果
│   ├── CP-*.csv               # CP 方法模拟结果
│   ├── CP2-*.csv              # CP2 方法模拟结果
│   ├── EL-*.csv               # EL 方法模拟结果
│   └── PCEL*.csv              # PCEL 方法模拟结果
└── eps/                        # 实验图片结果（EPS格式）
    ├── fig4-a.eps             # 图4-a
    ├── fig4-b.eps             # 图4-b
    └── PCEL*.eps              # PCEL 方法相关图片
```

## 依赖项

运行本代码需要安装以下 R 包：

```r
install.packages(c("spdep", "emplik", "expm", "maxLik", "sp", "terra", 
                   "MASS", "sf", "spData"))
```

## 主要文件说明

### RealSimPCEL.R

使用 PCEL 方法对 SARAR 模型进行参数估计和假设检验的代码。主要功能包括：

- **数据来源**: Columbus 犯罪数据（`columbus` 数据集，来自 `spdep` 包）
- **模型设定**: 
  - 样本量 `n = 49`
  - 回归系数维度 `p = 3`
  - 空间权重矩阵 `Wn` 和 `Mn`（基于 `col.gal.nb` 邻接关系）
  - 参数：`beta = c(0, -0.5, -0.5)`, `rou1 = 0.5`, `rou2 = 0.5`, `sigma2 = 1`
- **主要函数**:
  - `myPCEL()`: 计算 PCEL 统计量
  - `Yn_hat()`: 预测响应变量
  - `f()`: 矩阵向量乘积辅助函数
  - `tr()`: 矩阵迹计算函数
  - `L2sqare()`: 矩阵对角线元素平方和函数
- **输出**:
  - 参数估计值
  - PCEL 统计量
  - 假设检验结果（p 值）
  - 平均绝对误差（MAE）

### RealSimEL.R

使用传统 EL 方法对 SARAR 模型进行参数估计的代码，与 `RealSimPCEL.R` 类似但使用不同的估计方法。

### PCAEL.R

PCA 经验似然统计量计算脚本，用于生成模拟数据并进行统计分析。

### 核心函数（core/ 目录）

- **core/GlambdaChen.R**: 计算实证似然方法中的 Lagrange 乘数（Lambda）的核心函数
- **core/PCEL.R**: PCEL 方法的核心实现函数

### 数据处理和可视化脚本

- **PCA_Data.R / PCA_Data2.R**: PCA 数据处理脚本
- **PCA_latex.R / PCA_latex2.R**: 生成 LaTeX 格式的输出表格
- **PCA_QQ.R**: 绘制 Q-Q 图进行分布检验
- **fig4-ab.R**: 生成论文图4（a和b）
- **fig5-abc.R**: 生成论文图5（a、b和c）

## 使用方法

### 快速开始

1. **设置工作目录**：在 R 中设置工作目录到项目根目录
   
   ```r
   setwd("path/to/PCEL")
   ```
   
2. **运行 PCEL 方法**：
   ```r
   source("core/PCEL.R")
   ```


### 修改参数

在 `code/RealSimPCEL.R` 文件的第 41-55 行可以修改以下参数：

- `n`: 样本量（默认：49）
- `p`: 回归系数维度（默认：3）
- `beta`: 回归系数真值（默认：`c(0, -0.5, -0.5)`）
- `rou1`, `rou2`: 空间自回归系数（默认：0.5）
- `sigma2`: 误差方差（默认：1）
- `mu3`, `mu4`: 三阶和四阶矩
- `s`: PCEL 方法中的主成分数量（约第 129 行）

### 参数估计

代码使用 `nlminb()` 函数进行优化，参数约束为：
- `beta`: 无约束（第一个元素）或非正约束（后两个元素）
- `rou1`, `rou2`: [-0.99, 0.99]
- `sigma2`: [0.01, Inf]
- `mu3`, `mu4`: 无约束

### 运行模拟实验

运行 `code/PCAEL.R` 中的 `SARAR_Data()` 函数可以进行大规模模拟实验：

```r
source("code/PCAEL.R")
# 示例：m=10 (n=100), error=1, indexs=c(0,0.1), nsim=1000
results <- SARAR_Data(m=10, error=1, indexs=c(0,0.1), nsim=1000)
```

## 数据文件

`data/` 目录包含不同方法在不同样本量和场景下的模拟结果：

- **CP-*.csv**: CP 方法模拟结果（场景1-4，样本量：100, 225, 400, 625, 900）
- **CP2-*.csv**: CP2 方法模拟结果（场景1，样本量：100, 225, 400, 625, 900）
- **EL-*.csv**: 传统实证似然方法的结果（场景1-4，样本量：100, 225, 324, 400, 625, 900）
- **PCEL1-*.csv**: PCEL 方法结果（场景1，样本量：100, 225, 400, 625, 900）
- **PCEL2-*.csv**: PCEL 方法结果（场景2，样本量：100, 225, 400, 625, 900）

文件命名格式：`方法-场景-样本量.csv`

例如：`PCEL1-1-100.csv` 表示 PCEL 方法、场景1、样本量100的模拟结果。

## 注意事项

1. **路径设置**：确保 `GlambdaChen.R` 文件在正确路径下，或修改 `source()` 语句中的路径
   - `code/RealSimPCEL.R` 和 `code/RealSimEL.R` 使用 `source('GlambdaChen.R')`
   - `code/PCAEL.R` 和 `core/PCEL.R` 也依赖 `GlambdaChen.R`

2. **空间权重矩阵**：`Wn` 和 `Mn` 需要满足模型假设，通常基于空间邻接关系构建

3. **优化算法**：使用 `nlminb()` 进行优化，可能需要调整初始值以获得更好的收敛性

4. **R 项目文件**：可以使用 `code/PCEL.Rproj` 在 RStudio 中打开项目

5. **数据文件**：Columbus 空间数据文件（`.dbf`, `.shp`, `.shx`）位于 `code/` 目录下

6. **结果文件**：模拟结果保存在 `data/` 目录，图片结果保存在 `eps/` 目录

## 输出说明

运行 `code/RealSimPCEL.R` 后会输出：

1. **参数初值**: 优化算法的初始参数值
2. **参数估计值**: 优化后的参数估计
3. **假设检验结果**: PCEL 统计量是否小于卡方分布的临界值
4. **MAE**: 平均绝对误差
5. **p 值**: 假设检验的 p 值

运行 `code/PCAEL.R` 中的模拟函数会生成 CSV 格式的结果文件，保存在 `data/` 目录下。

运行可视化脚本（如 `fig4-ab.R`, `fig5-abc.R`）会生成 EPS 格式的图片文件，保存在 `eps/` 目录下。

## 许可证

本项目代码仅供学术研究使用。
