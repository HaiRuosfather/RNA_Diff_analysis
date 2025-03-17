# RNA-Seq差异分析脚本

## 简介
这是一个用于RNA-Seq数据差异分析的脚本工具，支持从数据准备到差异基因分析的全流程操作。由简少开发，版本1.0发布于2025年3月14日。

## 安装与配置

### 1. 环境配置
1. **安装conda**：请参考[conda官方文档](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)进行安装。
2. **创建运行环境**：
   ```bash
   bash 010.create_env.sh
   ```
3. **安装依赖包**：
   ```bash
   bash 011.install_pkg.sh
   ```
4. **调整脚本路径**：在`100.RNA_seq.R`中，确保`source`函数路径指向`jzx_Function`库中的函数。

### 2. 配置文件
配置文件`DEG_config.tsv`是运行分析的核心。请根据以下参数进行配置：

| 参数 | 说明 |
|------|------|
| GSE_number | 芯片数据的GSE号（可选） |
| species | 物种信息（支持人、小鼠、大鼠） |
| projeid | 项目编号 |
| plotid | 输出文件编号 |
| ID_type | 原始数据的基因ID类型（如`SYMBOL`） |
| pvalue_column | 差异分析中p值的类型 |
| pvalue_cutoff | p值筛选标准 |
| FC_cutoff | 折变值（Fold Change）筛选标准 |
| go_analysis | 是否进行GO/KEGG分析（`TRUE`或`FALSE`） |
| analysis_type | 分析方法 |
| [sample_group] | 样本信息（样本名）和分组信息（分组名） |

**注意**：`sample_group`部分的样本信息可以与表达矩阵列名不一致，但分组信息必须按实验组在前、对照组在后的顺序排列。

### 3. 数据准备
- 数据文件格式：以`tab`分隔的文本文件，文件名为`[项目编号].txt`或`[项目编号].tsv`。
- 数据结构：行名为基因ID，列名为样本名称。列名排序无要求，脚本会自动与配置文件中的样本名取交集。

## 运行分析

1. **将文件放置到输出路径**：
   将配置文件、数据文件放置到指定的输出路径。
2. **激活运行环境**：
   ```bash
   conda activate R44
   ```
3. **运行脚本**：
   ```bash
   Rscript scripts/100.RNA_seq.R
   ```

## 输出结果
运行完成后，结果将保存在当前目录下，具体文件格式和内容请参考脚本说明。

## 联系方式
如在使用过程中遇到问题，请联系开发者：[简仲翔](mailto:1416684627@qq.com)。

---

**注意**：请根据实际情况调整路径和联系方式等信息。
