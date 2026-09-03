
实现一个脚本，用来比对 smc 策略变化对于 其输出序列的 head 和 tail 的影响有多大。下面是详细说明

接收两个 fastq 文件。接收两组文件夹。一个是实验组，一个是对照组。
基于文件夹中的 文件名 将 两组实验的 对应结果关联起来。
文件是 fastq 文件，基于文件中的 qname 将文件中的 序列粒度的结果对应起来。
然后取 序列对 的 head 50 ,tail 50. 然后计算 edit-distance. 然后 head 的 editdistance 和 tail 的 edit distance 分别收集起来。后面用来做 柱状图。
然后也增加一些统计指标。可以让我判断，这次 策略变更是否可以上线

代码写到@/root/projects/gsda/gseda/src/gseda/smc_ab_test/head_tail_ab_analysis.py 中。