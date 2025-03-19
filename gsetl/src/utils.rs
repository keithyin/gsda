use rand::{self, Rng};

fn partition(data: &mut [f64], pivot_index: usize) -> usize {
    data.swap(pivot_index, data.len() - 1);
    let pivot = data[data.len() - 1];
    
    let mut store_index = 0;
    for i in 0..data.len() - 1 {
        if data[i] <= pivot {
            data.swap(i, store_index);
            store_index += 1;
        }
    }
    
    data.swap(store_index, data.len() - 1);
    store_index
}

fn select_nth(data: &mut [f64], n: usize) -> f64 {
    let mut left = 0;
    let mut right = data.len() - 1;
    let mut rng = rand::rng();
    
    loop {
        if left == right {
            return data[left];
        }
        
        let pivot_index = rng.random_range(left..=right);
        let pivot_index = partition(&mut data[left..=right], pivot_index - left) + left;

        if n == pivot_index {
            return data[n];
        } else if n < pivot_index {
            right = pivot_index - 1;
        } else {
            left = pivot_index + 1;
        }
    }
}

pub fn median(data: &mut Vec<f64>) -> Option<f64> {
    if data.is_empty() {
        return None;
    }
    
    let len = data.len();
    let mid = len / 2;
    
    if len % 2 == 0 {
        // 偶数个元素：找出中间两个数并返回平均值
        let left_mid = select_nth(data, mid - 1);
        let right_mid = select_nth(data, mid);
        Some((left_mid + right_mid) / 2.0)
    } else {
        // 奇数个元素：找出中间的数
        Some(select_nth(data, mid))
    }
}