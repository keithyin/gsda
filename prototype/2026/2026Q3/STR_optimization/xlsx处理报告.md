# 三批 xlsx（样本 ↔ barcode ↔ RUN）处理报告

- 统一表：`样本-run-barcode对应关系.tsv` / `样本-run-barcode对应关系.xlsx`（共 68 行）
- barcode 映射（24标签 ↔ 新编号 ↔ 原编号）来自 数据介绍.md 的映射表

## 第一批  (`first-batch-of-data/`)
- xlsx：`第一批-样本对应barcode.xlsx`，共 10 行样本记录
- 盘上 run 目录：20260720_250302Y0001_Run0003, 20260722_250302Y0006_Run0002, 20260723_250302Y0001_Run0003
- block 划分（RUN 号所在行起为新 block）：
  - （xlsx 未标注 RUN 号） | DPN=- | 服务器=- | 10 个样本
- 参考序列目录：`first-batch-of-data/STR第一批一代测序/STR一代测序/merged_output`（独立文件 10 / 汇总文件内 0 / 缺失 0）
- xlsx 无 RUN 号，经确认仅采用：20260723_250302Y0001_Run0003
### run: `20260723_250302Y0001_Run0003`

## 第二批  (`second-batch-of-data/`)
- xlsx：`第二批-样本对应barcode.xlsx`，共 44 行样本记录
- 盘上 run 目录：20260805_250302Y0001_Run0007, 20260805_250804Y0004_Run0001, 20260805_250804Y0004_Run0002
- block 划分（RUN 号所在行起为新 block）：
  - 20260805_250804Y0004_Run0001 | DPN=1-15 | 服务器=153服务器-PRO | 15 个样本
  - 20260805_250804Y0004_Run0002 | DPN=16-29 | 服务器=153服务器-PRO | 14 个样本
  - 20260805_250302Y0001_Run0007 | DPN=1-15 | 服务器=72 服务器 | 15 个样本
- 参考序列目录：`second-batch-of-data/sanger_result/merged_output`（独立文件 26 / 汇总文件内 3 / 缺失 0）
### run: `20260805_250302Y0001_Run0007`

### run: `20260805_250804Y0004_Run0001`
- 盘上未分配（池内位置未用于本 run 的 block）：24标签-16(Adaptor-barcode217-2, 2reads)
- 池外 barcode 有 reads（疑似邻码误分）：`Adaptor-barcode215-4`(1reads), `Adaptor-barcode222-4`(2reads), `Adaptor-barcode222-5`(3reads), `Adaptor-barcode230-4`(2reads), `Adaptor-barcode230-5`(1reads), `Adaptor-barcode251-2`(1reads), `Adaptor-barcode253-2`(2reads), `Adaptor-barcode254-2`(2reads), `Adaptor-barcode255-2`(1reads), `Adaptor-barcode256-5`(1reads), `Adaptor-barcode264-3`(1reads)

### run: `20260805_250804Y0004_Run0002`
- 盘上未分配（池内位置未用于本 run 的 block）：24标签-15(Adaptor-barcode253-4, 3reads), 24标签-16(Adaptor-barcode217-2, 2reads)
- 池外 barcode 有 reads（疑似邻码误分）：`Adaptor-barcode201-2`(1reads), `Adaptor-barcode222-4`(3reads), `Adaptor-barcode222-5`(1reads), `Adaptor-barcode230-4`(2reads), `Adaptor-barcode230-5`(2reads), `Adaptor-barcode250-2`(1reads), `Adaptor-barcode254-2`(1reads), `Adaptor-barcode266-4`(1reads), `Adaptor-barcode266-5`(2reads)

## 第三批  (`third-batch-of-data/`)
- xlsx：`第三批-样本对应barcode 30-43.xlsx`，共 14 行样本记录
- 盘上 run 目录：20260815_250804Y0004_Run0002
- block 划分（RUN 号所在行起为新 block）：
  - 20260815_250804Y0004_Run0002 | DPN=30-43 | 服务器=153服务器 | 14 个样本
- 参考序列目录：`third-batch-of-data/STR第三批一代测序/STR第三批一代测序/merged_output`（独立文件 14 / 汇总文件内 0 / 缺失 0）
- 有参考序列但 ccs xlsx 无对应样本：260805STR16-1.fa（16-1）
### run: `20260815_250804Y0004_Run0002`
- 盘上未分配（池内位置未用于本 run 的 block）：24标签-23(Adaptor-barcode233-0, 1reads)
- 池外 barcode 有 reads（疑似邻码误分）：`Adaptor-barcode201-0`(1reads), `Adaptor-barcode202-0`(4reads), `Adaptor-barcode204-0`(1reads), `Adaptor-barcode205-0`(1reads), `Adaptor-barcode210-0`(1reads), `Adaptor-barcode213-0`(3reads), `Adaptor-barcode214-0`(1reads), `Adaptor-barcode215-0`(1reads), `Adaptor-barcode221-0`(1reads), `Adaptor-barcode222-0`(1reads), `Adaptor-barcode222-4`(2reads), `Adaptor-barcode222-5`(6reads), `Adaptor-barcode225-0`(1reads), `Adaptor-barcode229-4`(1reads), `Adaptor-barcode230-4`(2reads), `Adaptor-barcode232-0`(4reads), `Adaptor-barcode234-0`(1reads), `Adaptor-barcode235-0`(1reads), `Adaptor-barcode236-0`(1reads), `Adaptor-barcode237-4`(1reads), `Adaptor-barcode239-0`(1reads), `Adaptor-barcode241-0`(1reads), `Adaptor-barcode244-0`(1reads), `Adaptor-barcode245-0`(1reads), `Adaptor-barcode246-0`(2reads), `Adaptor-barcode247-0`(3reads), `Adaptor-barcode248-0`(1reads), `Adaptor-barcode250-0`(2reads), `Adaptor-barcode251-0`(4reads), `Adaptor-barcode252-0`(2reads), `Adaptor-barcode253-0`(2reads), `Adaptor-barcode253-2`(3reads), `Adaptor-barcode254-0`(1reads), `Adaptor-barcode254-2`(2reads), `Adaptor-barcode255-0`(1reads), `Adaptor-barcode259-0`(1reads), `Adaptor-barcode260-0`(2reads), `Adaptor-barcode261-4`(1reads), `Adaptor-barcode262-0`(1reads), `Adaptor-barcode266-0`(1reads), `Adaptor-barcode266-4`(2reads), `Adaptor-barcode266-5`(3reads), `Adaptor-barcode267-0`(1reads), `Adaptor-barcode268-0`(2reads), `Adaptor-barcode270-4`(1reads), `Adaptor-barcode271-0`(1reads), `Adaptor-barcode272-0`(2reads), `Adaptor-barcode274-0`(2reads), `Adaptor-barcode279-0`(2reads), `Adaptor-barcode282-0`(1reads), `Adaptor-barcode283-0`(2reads), `Adaptor-barcode286-0`(1reads), `Adaptor-barcode286-3`(1reads), `Adaptor-barcode289-0`(1reads), `Adaptor-barcode289-4`(1reads), `Adaptor-barcode291-0`(2reads), `Adaptor-barcode292-0`(2reads), `Adaptor-barcode294-0`(1reads), `Adaptor-barcode295-0`(2reads)

