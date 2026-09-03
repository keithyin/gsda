smc 输出的bam包含 di 字段。
该字段是一个 string。但是实际包含的信息是一个list，使用 ";" 分隔。
每个单元由 pos,base,frac,depth,phreq 构成。
di 的每个元素表示的是 subreads 在该位点有支撑，但是 smc输出 seq 上并没有该 base。
需要注意，pos 可能重复，意味着这个位置存在多个 subreads 有支撑，但是 smc seq 没有体现的 base。
di中的 qual 是 -10*log(1-p), 没有其它处理 

写一个脚本，画出 di 中的 qual 和 smc seq qual 的分布，使用 histogram。然后使用 legend 区分 两组数据。
代码需要一些过滤逻辑，比如 rq 过滤。np 过滤。rq 读取 bam record rq 字段，np读取 bam record 的 np 字段。

输入文件，参考 @/data1/ccs_data/ecoli-data-for-marketing/20250623_250302Y0002_Run0004_adapter.withdel.smc_all_reads.bam
生成的 图是  *.di_vs_seq.png

请实现对应的代码，遵循 SOLID 原则，并不要过度抽象。并注意，输入、输出 bam 并不是 比对的 结果 bam。
代码写到  @/root/projects/gsda/gseda/src/gseda/smc_bam_post_porcess/di_vs_seq_qual.py 中。