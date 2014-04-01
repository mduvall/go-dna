package dna

type Dna struct {
	Sequence []byte
}

func codonTranslationMap() map[string]string {
	codonTranslation := map[string]string{
		"UUU": "F",
		"CUU": "L",
		"AUU": "I",
		"GUU": "V",
		"UUC": "F",
		"CUC": "L",
		"AUC": "I",
		"GUC": "V",
		"UUA": "L",
		"CUA": "L",
		"AUA": "I",
		"GUA": "V",
		"UUG": "L",
		"CUG": "L",
		"AUG": "M",
		"GUG": "V",
		"UCU": "S",
		"CCU": "P",
		"ACU": "T",
		"GCU": "A",
		"UCC": "S",
		"CCC": "P",
		"ACC": "T",
		"GCC": "A",
		"UCA": "S",
		"CCA": "P",
		"ACA": "T",
		"GCA": "A",
		"UCG": "S",
		"CCG": "P",
		"ACG": "T",
		"GCG": "A",
		"UAU": "Y",
		"CAU": "H",
		"AAU": "N",
		"GAU": "D",
		"UAC": "Y",
		"CAC": "H",
		"AAC": "N",
		"GAC": "D",
		"UAA": "Stop",
		"CAA": "Q",
		"AAA": "K",
		"GAA": "E",
		"UAG": "Stop",
		"CAG": "Q",
		"AAG": "K",
		"GAG": "E",
		"UGU": "C",
		"CGU": "R",
		"AGU": "S",
		"GGU": "G",
		"UGC": "C",
		"CGC": "R",
		"AGC": "S",
		"GGC": "G",
		"UGA": "Stop",
		"CGA": "R",
		"AGA": "R",
		"GGA": "G",
		"UGG": "W",
		"CGG": "R",
		"AGG": "R",
		"GGG": "G",
	}

	return codonTranslation
}

func nucleobaseComplementMap() map[byte]byte {
	return map[byte]byte{
		'A': 'T',
		'T': 'A',
		'C': 'G',
		'G': 'C',
	}
}

func (dna *Dna) complement(reverse bool) []byte {
	translationTable := nucleobaseComplementMap()
	translation := make([]byte, len(dna.Sequence))

	if reverse {
		for i, val := range dna.Sequence {
			translation[len(dna.Sequence)-i-1] = translationTable[val]
		}
	} else {
		for i, val := range dna.Sequence {
			translation[i] = translationTable[val]
		}
	}

	return translation
}

func (dna *Dna) ReverseComplement() []byte {
	return dna.complement(true)
}

func (dna *Dna) Complement() []byte {
	return dna.complement(false)
}

func (dna *Dna) GetRnaTranslationSequence() []byte {
	translation := make([]byte, len(dna.Sequence))

	for i, val := range dna.Sequence {
		if val == 'T' {
			translation[i] = 'U'
		} else {
			translation[i] = val
		}
	}

	return translation
}

func (dna *Dna) GetCodonSequence() []string {
	codonSequence := make([]string, 0)
	rnaSequence := dna.GetRnaTranslationSequence()

	for i := 0; i < len(rnaSequence)-(len(rnaSequence)%3); i += 3 {
		codonSequence = append(codonSequence, string(rnaSequence[i:i+3]))
	}

	return codonSequence
}

func (dna *Dna) GetAminoAcidSequence() []string {
	codonTranslationMap := codonTranslationMap()
	aminoAcidSequence := make([]string, 0)

	for _, codon := range dna.GetCodonSequence() {
		aminoAcidSequence = append(aminoAcidSequence, codonTranslationMap[codon])
	}

	return aminoAcidSequence
}

func (dna *Dna) GetOpenReadingFrames() []Dna {
	orf := make([]Dna, 6)

	// The three reading frames are shifted by 1 nucleobase up until another
	// codon is valid, which is 3 positions later (same as [0:] sequence)
	orf[0] = Dna{dna.Sequence[0:]}
	orf[1] = Dna{dna.Sequence[1:]}
	orf[2] = Dna{dna.Sequence[2:]}

	// The next three reading frames are the same as above except the reverse
	// complement of the DNA strand
	reverseComplementSequence := dna.ReverseComplement()
	orf[3] = Dna{reverseComplementSequence[0:]}
	orf[4] = Dna{reverseComplementSequence[1:]}
	orf[5] = Dna{reverseComplementSequence[2:]}

	return orf
}
