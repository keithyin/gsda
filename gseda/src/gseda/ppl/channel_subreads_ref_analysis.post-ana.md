一个文件夹，里面有很多子目录 以 -filtered-metric 结尾。
filtered-metric文件夹中有一个以 .gsmm2-hp-aggr.csv 结尾的文件

字段含义如下：
motif: (A)3 表示 AAA, (A)4 表示 AAAA, 其余类似。motif 表示是 reference
true_base: 表示 reference 的 base
true_cnt: 表示 poly-N 的重复次数， motif 就是用 (true_base)true_cnt构成的
called: 该 poly-N 区域 call 出来多少个 base
tag: mixed 表示 该 poly-N 区域不仅有 和 true_base 一致的base，还存在不一致的 base

现在我想要做一下数据分析： 使用 polars
1. 将所有 gsmm2-hp-aggr.csv 的 tag == pure+mixed 的数据过滤出来
2. 对于每个 motif，仅考虑 called in [true_cnt-1, true_cnt, true_cnt]。
3. 开始绘图。不同的 true_base 放到不同的 .jpg 图像中，同一个 true_base ，相同的 motif 的放在一axe 中，其中不同的 called 作为不同的 legend。然后对于 ratio_within_motif_tag 画柱状图。
4. 所有的绘图结果保存。