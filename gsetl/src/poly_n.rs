pub fn find_homopolymer_regions(sequence: &[u8]) -> Vec<(usize, usize, u8)> {
    if sequence.is_empty() {
        return Vec::new();
    }

    let mut regions = Vec::new();
    let mut start = 0;
    let mut bases = sequence.iter();
    let mut pre_base = *bases.next().unwrap();

    for (mut i, &cur_base) in bases.enumerate() {
        i += 1;
        if cur_base != pre_base {
            if (i - start) > 1 {
                regions.push((start, i, pre_base));
            }
            start = i;
            pre_base = cur_base;
        }
    }

    // 保存最后一个 homopolymer 区域
    if sequence.len() - start > 1 {
        regions.push((start, sequence.len(), pre_base));
    }

    regions
}

#[derive(Debug, Clone, Copy)]
pub enum PosRelation {
    Left,
    Right,
    Middle,
}

pub fn position_relation(seb: &(usize, usize, u8), pos: usize) -> PosRelation {
    return if pos < seb.0 {
        PosRelation::Left
    } else if pos >= seb.1 {
        PosRelation::Right
    } else {
        PosRelation::Middle
    };
}

#[cfg(test)]
mod test {
    use crate::poly_n::find_homopolymer_regions;

    #[test]
    fn test_find_homopolymer_regions() {
        let seq = b"AAAACCCGGTT";
        let res = find_homopolymer_regions(seq);
        assert_eq!(res, vec![(0, 4, 65), (4, 7, 67), (7, 9, 71), (9, 11, 84)]);

        let seq = b"AAAACCCGGT";
        let res = find_homopolymer_regions(seq);
        assert_eq!(res, vec![(0, 4, 65), (4, 7, 67), (7, 9, 71)]);

        let seq = b"ACCCGGT";
        let res = find_homopolymer_regions(seq);
        assert_eq!(res, vec![(1, 4, 67), (4, 6, 71)]);

        let seq = b"AACGGT";
        let res = find_homopolymer_regions(seq);
        // println!("{:?}", res);
        assert_eq!(res, vec![(0, 2, 67), (3, 5, 71)]);
    }
}
