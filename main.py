class DNAStorage:
    def __init__(self, max_homopolymer=2):

        # Define transition rules to prevent homopolymers and control nucleotide usage
        # When the last two nucleotides are the same, use these rules
        self.transition_rule_double = {
            'A': ['C', 'G', 'T'],  # Exclude 'A' to prevent 'AAA'
            'C': ['A', 'T'],       # Exclude 'C' to prevent 'CCC'
            'G': ['A', 'T'],       # Exclude 'G' to prevent 'GGG'
            'T': ['A', 'C', 'G'],  # Exclude 'T' to prevent 'TTT'
        }

        # When the last two nucleotides are not the same, use these rules
        self.transition_rules = {
            'A': ['A', 'C', 'G', 'T'],    # All nucleotides can follow 'A'
            'C': ['A', 'T'],              # Reduce 'G' and 'C' after 'C'
            'G': ['A', 'T'],              # Reduce 'G' and 'C' after 'G'
            'T': ['A', 'C', 'G', 'T'],    # All nucleotides can follow 'T'
            '': ['A', 'C', 'G', 'T']      # For the first nucleotide
        }

    def string_to_binary(self, text):
        binary = ''.join(format(ord(char), '08b') for char in text)
        return binary

    def binary_to_string(self, binary):
        # If the binary string length isn't a multiple of 8, pad it with zeros on the right
        if len(binary) % 8 != 0:
            binary = binary.ljust(len(binary) + (8 - len(binary) % 8), '0')

        text = ''
        for i in range(0, len(binary), 8):
            byte = binary[i:i + 8]     # Extract 8 bits
            text += chr(int(byte, 2))  # Convert to character and append
        return text

    def binary_to_dna(self, binary):

        dna = ''
        prev_nucleotide = ''
        i = 0
        # Process the binary string until all bits are converted
        while i < len(binary):
            # Determine which transition rules to use based on the last two nucleotides
            if len(dna) >= 2 and dna[-1] == dna[-2]:
                # Last two nucleotides are the same; use transition_rule_double
                possible_nucleotides = self.transition_rule_double[prev_nucleotide]
            else:
                # Use regular transition_rules
                possible_nucleotides = self.transition_rules[prev_nucleotide]

            # Get the number of possible nucleotides for the current position (could be 2 or 4)
            num_nucleotides = len(possible_nucleotides)

            # Determine the number of bits to process based on available nucleotides
            if num_nucleotides == 4:
                bits_to_process = 2
            elif num_nucleotides >= 2:
                bits_to_process = 1

            bits = binary[i:i + bits_to_process]
            # If there aren't enough bits left, pad with zeros on the right
            if len(bits) < bits_to_process:
                bits = bits.ljust(bits_to_process, '0')

            index = int(bits, 2)

            # Ensure the index is within the range of possible nucleotides
            if index >= num_nucleotides:
                index = index % num_nucleotides  # Wrap around using modulo

            # Select the nucleotide corresponding to the index
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

            # Convert the index back to a binary format string with leading zeros
            bits = format(index, f'0{bits_to_process}b')
            # Append the bits to the binary string
            binary += bits
            prev_nucleotide = nucleotide
            i += 1
        return binary

    def encode(self, text):
        binary = self.string_to_binary(text)
        dna = self.binary_to_dna(binary)
        return dna

    def decode(self, dna_sequence):
        try:
            binary = self.dna_to_binary(dna_sequence)
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

    encoded_dna = dna_storage.encode(original_text)
    print(f"Encoded DNA:\n{encoded_dna}\n")

    dna_storage.analyze_dna(encoded_dna)
    print()

    decoded_text = dna_storage.decode(encoded_dna)
    print(f"Decoded text: {decoded_text}")

if __name__ == "__main__":
    main()
