smc 输出的bam包含 di 字段。
该字段是一个 string。但是实际包含的信息是一个list，使用 ";" 分隔。
每个单元由 pos,base,frac,depth,phreq 构成。
di 的每个元素表示的是 subreads 在该位点有支撑，但是 smc输出 seq 上并没有该 base。
需要注意，pos 可能重复，意味着这个位置存在多个 subreads 有支撑，但是 smc seq 没有体现的 base。
现在写一个代码，将 di 中的 base 插入到 smc seq 中。其中 base 使用 di 中的 base。qual 使用 di 中的 phreq（di中的直接是 -10*log(1-p), 没有其它处理）
插入之后，重新计算整个 reads 的 q。然后更新 rq 字段的值。

输入文件，参考 @/data1/ccs_data/ecoli-data-for-marketing/20250623_250302Y0002_Run0004_adapter.withdel.smc_all_reads.bam
生成的文件名为 *.with_di.bam

请实现对应的代码，遵循 SOLID 原则，并不要过度抽象。并注意，输入、输出 bam 并不是 比对的 结果 bam。
代码写到  @/root/projects/gsda/gseda/src/gseda/smc_bam_post_porcess/insert_di.py 中。