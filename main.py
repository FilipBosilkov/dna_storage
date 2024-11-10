class DNAStorage:
    def __init__(self, max_homopolymer=2):

        self.transition_rule_double = {
            'A': ['C', 'G', 'T'],
            'C': ['A', 'T'],
            'G': ['A', 'T'],
            'T': ['A', 'C', 'G'],
        }

        self.transition_rules = {
            'A': ['A', 'C', 'G', 'T'],
            'C': ['A', 'T'],
            'G': ['A', 'T'],
            'T': ['A', 'C', 'G', 'T'],
            '': ['A', 'C', 'G', 'T']
        }

    def string_to_binary(self, text):
        binary = ''.join(format(ord(char), '08b') for char in text)
        return binary

    def binary_to_string(self, binary):
        if len(binary) % 8 != 0:
            binary = binary.ljust(len(binary) + (8 - len(binary) % 8), '0')

        text = ''
        for i in range(0, len(binary), 8):
            byte = binary[i:i + 8]
            text += chr(int(byte, 2))
        return text

    def binary_to_dna(self, binary):

        dna = ''
        prev_nucleotide = ''
        i = 0
        while i < len(binary):
            if len(dna) >= 2 and dna[-1] == dna[-2]:
                possible_nucleotides = self.transition_rule_double[prev_nucleotide]
            else:
                possible_nucleotides = self.transition_rules[prev_nucleotide]

            num_nucleotides = len(possible_nucleotides)

            if num_nucleotides == 4:
                bits_to_process = 2
            elif num_nucleotides >= 2:
                bits_to_process = 1

            bits = binary[i:i + bits_to_process]
            if len(bits) < bits_to_process:
                bits = bits.ljust(bits_to_process, '0')

            index = int(bits, 2)

            if index >= num_nucleotides:
                index = index % num_nucleotides

            nucleotide = possible_nucleotides[index]
            dna += nucleotide
            prev_nucleotide = nucleotide
            i += bits_to_process
        return dna

    def dna_to_binary(self, dna):

        binary = ''
        prev_nucleotide = ''
        i = 0
        while i < len(dna):
            nucleotide = dna[i]
            if i >= 2 and dna[i - 1] == dna[i - 2]:
                possible_nucleotides = self.transition_rule_double[prev_nucleotide]
            else:
                possible_nucleotides = self.transition_rules[prev_nucleotide]

            num_nucleotides = len(possible_nucleotides)

            if num_nucleotides == 4:
                bits_to_process = 2
            elif num_nucleotides >= 2:
                bits_to_process = 1

            index = possible_nucleotides.index(nucleotide)

            bits = format(index, f'0{bits_to_process}b')
            binary += bits
            prev_nucleotide = nucleotide
            i += 1
        return binary

    def encode(self, text, copies=4):
        binary = self.string_to_binary(text)
        dna = self.binary_to_dna(binary)
        prefix_separator = 'ATGC'
        suffix_separator = 'GCTA'
        final_sequence = (prefix_separator + dna + suffix_separator) * copies
        return final_sequence

    def decode(self, dna_sequence):
        copies = []
        sequences = dna_sequence.split('ATGC')
        for seq in sequences:
            if seq:
                seq = seq.split('GCTA')[0]
                if seq:
                    copies.append(seq)
        if not copies:
            return ""
        most_common = max(set(copies), key=copies.count)
        # most common sequence in case of loss
        try:
            binary = self.dna_to_binary(most_common)
            text = self.binary_to_string(binary)
            return text
        except ValueError as e:
            print(f"Decoding error: {e}")
            return ""

    def analyze_dna(self, dna_sequence):

        count_A = dna_sequence.count('A')
        count_C = dna_sequence.count('C')
        count_G = dna_sequence.count('G')
        count_T = dna_sequence.count('T')
        total_length = len(dna_sequence)

        if total_length == 0:
            print("The DNA sequence is empty.")
            return

        percentage_A = (count_A / total_length) * 100
        percentage_C = (count_C / total_length) * 100
        percentage_G = (count_G / total_length) * 100
        percentage_T = (count_T / total_length) * 100
        gc_content = (count_G + count_C) / total_length * 100

        print(f"Total nucleotides: {total_length}")
        print(f"Number of 'A's: {count_A} ({percentage_A:.2f}%)")
        print(f"Number of 'C's: {count_C} ({percentage_C:.2f}%)")
        print(f"Number of 'G's: {count_G} ({percentage_G:.2f}%)")
        print(f"Number of 'T's: {count_T} ({percentage_T:.2f}%)")
        print(f"GC Content: {gc_content:.2f}%")

def main():
    dna_storage = DNAStorage()

    original_text = "Hello! This is a test. Thank you."
    print(f"Original text: {original_text}\n")

    encoded_dna = dna_storage.encode(original_text, copies=4)
    print(f"Encoded DNA:\n{encoded_dna}\n")

    dna_storage.analyze_dna(encoded_dna)
    print()

    decoded_text = dna_storage.decode(encoded_dna)
    print(f"Decoded text: {decoded_text}")

if __name__ == "__main__":
    main()
