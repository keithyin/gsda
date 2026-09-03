

def main():
    seq1 = "CCCAAT-CGCT--CAG-TCAG--GT-AG-TGATG--CTTCT-TG"
    seq2 = "CCCAAT-CGCT--CAG-TCAG--GT-AG-TG-TG--C-----TG"
    
    count = 0
    eq = 0
    for i in range(len(seq1)):
        if seq1[i] == "-" and seq2[i] == "-":
            continue
        count += 1
        if seq1[i] == seq2[i]:
            eq += 1
    
    print(f"identity: {eq / count}")

if __name__ == "__main__":
    main()