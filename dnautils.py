from constants import *
from collections import Counter
import random

class DnaUtils:

    def __init__(self, seq, seq_type="DNA"):
        self.seq = seq
        self.seq_type = seq_type
        self.is_valid = self.validate()
        assert self.is_valid, f"The data is not valid for seq {self.seq}"

    def get_seq_info(self):
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Lenght]: {len(self.seq)}"

    def get_seq(self):
        return self.seq

    def validate(self):
        tmp_seq = self.seq.upper()
        for nuc in tmp_seq:
            if nuc not in NUCLEOTIDES[self.seq_type]:
                return False
        return tmp_seq
        #return set(NUCLEOTIDES[self.seq_type]).issuperset(self.seq)

    def nucleotide_frequemcy(self):
        freq_dic = {"A":0, "T":0, "C":0, "G":0, "U":0}
        for nuc in self.seq:
            freq_dic[nuc] += 1
        return freq_dic
    
    def transcription(self):
        if self.seq_type == "DNA":
            return self.seq.replace("T","U")
        return "Not DNA seq to transcript"
    
    def reverse_complement(self):
        """Used to get the other branch"""
        mapping = str.maketrans("ATCG","TAGC")
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        return round((self.seq.count('C') + self.seq.count('G'))/len(self.seq) * 100, 5)
    
    def gc_content_for_subseq(self, k=10):
        res =[]
        for i in range(0, len(self.seq)-k,k):
            subseq = self.seq[i:i+k]
            res.append(round((subseq.count('C') + subseq.count('G'))/len(subseq) * 100, 5))
        return res

    def codon_usage(self, aminoacid):
        tmpList =[]
        for i in range(0, len(self.seq) -2, 3):
            if self.seq == "DNA" and DNA_Codons[i] == aminoacid:
                tmpList.append(self.seq[i:i+3])
            if self.seq == "RNA" and RNA_Codons[i] == aminoacid:
                tmpList.append(self.seq[i:i+3])    
        freqDic = dict(Counter(tmpList))
        totalWight = sum(freqDic.values())
        for seq in freqDic:
            freqDic[seq] = round(freqDic[seq]/totalWight,2)
        return freqDic  


    def translate_to_aminoacids(self, initial_pos=0, reversed=False):
        amino_acids = []
        for pos in range(initial_pos, len(self.seq)-2, 3):
            current_seq = self.seq[pos:pos+3]
            if self.seq_type == "DNA":
                aa = DNA_Codons[current_seq]
            elif self.seq_type == "RNA":
                aa = RNA_Codons[current_seq]
            amino_acids.append(aa)
        if reversed: 
            return amino_acids[::-1]
        return amino_acids


    def gen_reading_frames(self):
        frames = []
        frames.append(self.translate_to_aminoacids(0))
        frames.append(self.translate_to_aminoacids(1))
        frames.append(self.translate_to_aminoacids(2))
        tmp_DnaUtils = DnaUtils(self.reverse_complement(), self.seq_type)
        frames.append(tmp_DnaUtils.translate_to_aminoacids(0))
        frames.append(tmp_DnaUtils.translate_to_aminoacids(1))
        frames.append(tmp_DnaUtils.translate_to_aminoacids(2))
        return frames
    
    def proteins_from_frame(self, frame):
        """frame = aminoacid sequence"""
        current_prot = []
        proteins = []
        for aa in frame:
            if aa == "_":#end_codon
                for p in current_prot:
                    proteins.append(p)
            if aa == "M":#start codon
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
        return proteins
    
    def count_kmer(self, kmer):
        kmer_count = 0
        for i in range(0, len(self.seq) - len(kmer) - 1):
            subseq = self.seq[i:len(kmer)]
            if subseq == kmer:
                kmer_count += 1
        return kmer_count
    
    def counter_kmer2(self, kmer):
        kmer_list = [self.seq[i:len(kmer)] for i in range(0, len(self.seq) - len(kmer) - 1) if self.seq[i:len(kmer)] == kmer]
        return len(kmer_list)
    
    def get_hamilton_distance_loop(self,seq2):
        hd = 0
        for i in range(len(self.seq)):
            if self.seq[i] != seq2[i]:
                hd += 1
        return hd

    def get_hamilton_distance_zip(self,seq2):
        zipped_dna = zip(self.seq,seq2)
        hd = [(n1, n2) for n1, n2 in zipped_dna if n1 != n2]
        return len(hd)
    
    