class RangeIdentityCalculator:
    def __init__(self, cigar_tuples):
        self.query_start_pos_and_op = []
        q_start_cursor = 0

        for op, length in cigar_tuples:
            self.query_start_pos_and_op.append((q_start_cursor, (op, length)))
            if op in {0, 1, 4, 7, 8}:  # Match, Ins, SoftClip, Equal, Diff
                q_start_cursor += length

    def range_identity(self, start, end):
        def binary_search(arr, key):
            low, high = 0, len(arr) - 1
            while low <= high:
                mid = (low + high) // 2
                if arr[mid][0] < key:
                    low = mid + 1
                elif arr[mid][0] > key:
                    high = mid - 1
                else:
                    return mid, True
            return low - 1, False

        start_idx, _ = binary_search(self.query_start_pos_and_op, start)
        end_idx, _ = binary_search(self.query_start_pos_and_op, end)

        # Ensure indices are within valid range
        start_idx = max(0, start_idx)
        end_idx = max(0, end_idx)

        res_cigars = []

        if start_idx == end_idx:
            op, length = self.query_start_pos_and_op[start_idx][1]
            start_pos = self.query_start_pos_and_op[start_idx][0]
            span_length = end - start
            overlap_start = max(start, start_pos)
            overlap_end = min(end, start_pos + length)

            if overlap_start >= overlap_end:
                return 0.0  # No valid overlap

            if op in {1, 8}:  # Insertion, Diff
                res_cigars.append((8, overlap_end - overlap_start))  # Treat as Diff
            elif op == 7:  # Equal
                res_cigars.append((7, overlap_end - overlap_start))
        else:
            # Handle the first segment
            start_pos = self.query_start_pos_and_op[start_idx][0]
            op, length = self.query_start_pos_and_op[start_idx][1]
            overlap_start = max(start, start_pos)
            overlap_end = min(end, start_pos + length)
            span_length = overlap_end - overlap_start

            if span_length > 0:
                if op in {1, 8}:  # Insertion, Diff
                    res_cigars.append((8, span_length))  # Treat as Diff
                elif op == 7:  # Equal
                    res_cigars.append((7, span_length))

            # Middle segments
            res_cigars.extend(
                (op, length)
                for _, (op, length) in self.query_start_pos_and_op[start_idx + 1 : end_idx]
            )

            # Handle the last segment
            end_pos = self.query_start_pos_and_op[end_idx][0]
            op, length = self.query_start_pos_and_op[end_idx][1]
            overlap_start = max(start, end_pos)
            overlap_end = min(end, end_pos + length)
            span_length = overlap_end - overlap_start

            if span_length > 0:
                if op in {1, 8}:  # Insertion, Diff
                    res_cigars.append((8, span_length))  # Treat as Diff
                elif op == 7:  # Equal
                    res_cigars.append((7, span_length))

        # Calculate spans and matches
        span_len = sum(length for op, length in res_cigars if op in {1, 7, 8})
        eq_len = sum(length for op, length in res_cigars if op == 7)

        return eq_len / (span_len if span_len > 0 else 1.0)
    
if __name__ == "__main__":
    # Example Usage
    # Example cigar string: [(0, 5), (8, 3), (7, 7), (1, 2)]
    # Operations: MATCH=0, INS=1, DIFF=8, EQUAL=7, etc.
    cigar_tuples = [
        (7, 5),  # eq
        (8, 3), # diff
        (7, 10), 
        (1, 2),
        (7, 6),
        (8, 4)
    ]

    calculator = RangeIdentityCalculator(cigar_tuples)
    identity = calculator.range_identity(0, 5)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(5, 8)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(8, 18)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(3, 7)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(7, 12)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(18, 26)
    print(f"Range Identity: {identity:.2f}")

    cigar_tuples = [
        (7, 5),
        (8, 5)
    ]

    calculator = RangeIdentityCalculator(cigar_tuples)
    identity = calculator.range_identity(0, 3)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(6, 10)
    print(f"Range Identity: {identity:.2f}")

    identity = calculator.range_identity(4, 7)
    print(f"Range Identity: {identity:.2f}")