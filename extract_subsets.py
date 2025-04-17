import numpy as np

def extract_fixed_length_rows(values: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """
    从每行中根据 start 和 end 提取固定长度的子序列。
    要求每行的 end - start 是相同的长度。
    """
    assert values.ndim == 2, "values 必须是二维数组"
    assert starts.shape == ends.shape and starts.ndim == 1, "starts 和 ends 应为相同长度的 1D 数组"
    
    lengths = ends - starts
    assert np.all(lengths == lengths[0]), "每行的提取长度必须一致"
    length = lengths[0]

    row_idx = np.arange(values.shape[0])[:, None]      # shape: [n, 1]
    col_idx = starts[:, None] + np.arange(length)      # shape: [n, length]

    return values[row_idx, col_idx]

if __name__ == "__main__":
    values = np.array([
        [1, 2, 3, 4, 5],
        [10, 20, 30, 40, 50],
        [100, 200, 300, 400, 500]
    ])  # shape = [3, 5]

    starts = np.array([1, 0, 2])
    ends = np.array([3, 2, 4])  # 每行都提取长度为 2 的切片

    # 调用函数
    result = extract_fixed_length_rows(values, starts, ends)

    # 输出结果
    print("提取结果：")
    print(result)