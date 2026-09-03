def rna_to_dna(rna_sequence):
    """
    将 RNA 序列转换为 DNA 序列
    原理：将所有的 'U' 替换为 'T'
    """
    return rna_sequence.upper().replace('U', 'T')

def main():
    print("RNA 到 DNA 序列转换器")
    print("=" * 30)
    
    # 获取用户输入
    rna_input = input("请输入 RNA 序列: ").strip()
    
    if not rna_input:
        print("错误：请输入有效的 RNA 序列")
        return
    
    # 转换为 DNA
    dna_sequence = rna_to_dna(rna_input)
    
    # 显示结果
    print(f"\n原始 RNA 序列: {rna_input.upper()}")
    print(f"转换后 DNA 序列: {dna_sequence}")
    
    # 验证转换
    u_count = rna_input.upper().count('U')
    t_count = dna_sequence.count('T')
    print(f"\n转换统计: 替换了 {u_count} 个 U 碱基")

if __name__ == "__main__":
    main()