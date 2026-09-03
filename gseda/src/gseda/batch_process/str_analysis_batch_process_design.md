
以下是 barcode 名字和新编号的对应关系
Adaptor-barcode251-4	barcode-01
Adaptor-barcode250-4	barcode-02
Adaptor-barcode271-3	barcode-03
Adaptor-barcode266-2	barcode-04
Adaptor-barcode229-0	barcode-05
Adaptor-barcode255-3	barcode-06
Adaptor-barcode261-3	barcode-07
Adaptor-barcode206-3	barcode-08
Adaptor-barcode230-2	barcode-09
Adaptor-barcode256-3	barcode-10
Adaptor-barcode207-4	barcode-11
Adaptor-barcode222-2	barcode-12
Adaptor-barcode270-0	barcode-13
Adaptor-barcode201-5	barcode-14
Adaptor-barcode253-4	barcode-15
Adaptor-barcode217-2	barcode-16
Adaptor-barcode214-1	barcode-17
Adaptor-barcode221-3	barcode-18
Adaptor-barcode216-0	barcode-19
Adaptor-barcode254-5	barcode-20
Adaptor-barcode220-0	barcode-21
Adaptor-barcode225-1	barcode-22
Adaptor-barcode233-0	barcode-23
Adaptor-barcode239-2	barcode-24


以下是(样本编号, barcode新编号别名, run号) 其实是标识了 (run号， barcode新编号别名) -> 样本编号的 对应关系。 barcode新编号别名和 barcode新编号的关系就是：  24标签-1 对应了 barcode01， 24标签-2 对应了 barcode02 ，以此类推
3N-1	24标签-1	20260805_250804Y0004_Run0001
5-3	24标签-2	20260805_250804Y0004_Run0001
5-8	24标签-3	20260805_250804Y0004_Run0001
5-10	24标签-4	20260805_250804Y0004_Run0001
5-12	24标签-5	20260805_250804Y0004_Run0001
5-15	24标签-6	20260805_250804Y0004_Run0001
7-10	24标签-7	20260805_250804Y0004_Run0001
12-9	24标签-8		20260805_250804Y0004_Run0001
25-1	24标签-9		20260805_250804Y0004_Run0001
27-1	24标签-10		20260805_250804Y0004_Run0001
27-4	24标签-11		20260805_250804Y0004_Run0001
27-9	24标签-12		20260805_250804Y0004_Run0001
28-1	24标签-13		20260805_250804Y0004_Run0001
28-2	24标签-14		20260805_250804Y0004_Run0001
28-3	24标签-15		20260805_250804Y0004_Run0001
29-1	24标签-1	20260805_250804Y0004_Run0002
29-2	24标签-2		20260805_250804Y0004_Run0002
29-3	24标签-3		20260805_250804Y0004_Run0002
29-5	24标签-4		20260805_250804Y0004_Run0002
29-8	24标签-5		20260805_250804Y0004_Run0002
29-9	24标签-6		20260805_250804Y0004_Run0002
29-10	24标签-7		20260805_250804Y0004_Run0002
33-3	24标签-8		20260805_250804Y0004_Run0002
33-5	24标签-9		20260805_250804Y0004_Run0002
33-6	24标签-10		20260805_250804Y0004_Run0002
33-11	24标签-11		20260805_250804Y0004_Run0002
38-1	24标签-12		20260805_250804Y0004_Run0002
38-2	24标签-13		20260805_250804Y0004_Run0002
38-3	24标签-14		20260805_250804Y0004_Run0002
3N-1	24标签-1	20260805_250302Y0001_Run0007
5-3	24标签-2		20260805_250302Y0001_Run0007
5-8	24标签-3		20260805_250302Y0001_Run0007
5-10	24标签-4		20260805_250302Y0001_Run0007
5-12	24标签-5		20260805_250302Y0001_Run0007
5-15	24标签-6		20260805_250302Y0001_Run0007
7-10	24标签-7		20260805_250302Y0001_Run0007
12-9	24标签-8		20260805_250302Y0001_Run0007
25-1	24标签-9		20260805_250302Y0001_Run0007
27-1	24标签-10		20260805_250302Y0001_Run0007
27-4	24标签-11		20260805_250302Y0001_Run0007
27-9	24标签-12		20260805_250302Y0001_Run0007
28-1	24标签-13		20260805_250302Y0001_Run0007
28-2	24标签-14		20260805_250302Y0001_Run0007
28-3	24标签-15		20260805_250302Y0001_Run0007


样本编号对应的真实 reference 序列在 @/data1/ccs_data/str-optimization/second_batch_of_data/sanger_result 目录下
测序结果在 @/data1/ccs_data/str-optimization/second_batch_of_data/20260805_*/barcode_assign/* 目录下

你先搞清楚上面的 测序结果 和 真实reference序列的对应关系，然后 调用
python /root/projects/gsda/gseda/src/gseda/ppl/reads_quality_stats_hp_v3.py --bams 所有fastq    --refs 对应参考序列  -f --ref-anchored --rq-range 0.99:1.1