module Rosalind

import Base.length

# (✓) Create abstract type for genetic sequence
# (✓) Add concrete protein sequence type
# ( ) Finish working on toProteinString

const DNAChars = ['A', 'C', 'G', 'T']
const RNAChars = ['A', 'C', 'G', 'U']
const ProteinChars = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
const RNACodonTable = ["UUU" => "F", "CUU" => "L", "AUU" => "I",
                       "GUU" => "V", "UUC" => "F", "CUC" => "L",
                       "AUC" => "I", "GUC" => "V", "UUA" => "L",
                       "CUA" => "L", "AUA" => "I", "GUA" => "V",
                       "UUG" => "L", "CUG" => "L", "AUG" => "M",
                       "GUG" => "V", "UCU" => "S", "CCU" => "P",
                       "ACU" => "T", "GCU" => "A", "UCC" => "S",
                       "CCC" => "P", "ACC" => "T", "GCC" => "A",
                       "UCA" => "S", "CCA" => "P", "ACA" => "T",
                       "GCA" => "A", "UCG" => "S", "CCG" => "P",
                       "ACG" => "T", "GCG" => "A", "UAU" => "Y",
                       "CAU" => "H", "AAU" => "N", "GAU" => "D",
                       "UAC" => "Y", "CAC" => "H", "AAC" => "N",
                       "GAC" => "D", "UAA" => "X", "CAA" => "Q",
                       "AAA" => "K", "GAA" => "E", "UAG" => "X",
                       "CAG" => "Q", "AAG" => "K", "GAG" => "E",
                       "UGU" => "C", "CGU" => "R", "AGU" => "S",
                       "GGU" => "G", "UGC" => "C", "CGC" => "R",
                       "AGC" => "S", "GGC" => "G", "UGA" => "X",
                       "CGA" => "R", "AGA" => "R", "GGA" => "G",
                       "UGG" => "W", "CGG" => "R", "AGG" => "R",
                       "GGG" => "G"]

abstract GeneticString

immutable DNAString <: GeneticString
  seq::String

  bases = Set(DNAChars)
  DNAString(seq) = fromAlphabet(bases, seq) ? new(seq) : error("String contains non-ACTG character.")
end

immutable RNAString <: GeneticString
  seq::String

  bases = Set(RNAChars)
  RNAString(seq) = fromAlphabet(bases, seq) ? new(seq) : error("String contains non-ACGU character.")
end

immutable ProteinString <: GeneticString
  seq::String

  aminos = Set(ProteinChars)
  ProteinString(seq) = fromAlphabet(aminos, seq) ? new(seq) : error("String contains invalid amino acid character")
end

function fromAlphabet(alphabet::Set{Char}, seq::String)
  all(b -> b in alphabet, seq)
end

function alphabet(dna::DNAString)
  DNAChars
end

function alphabet(rna::RNAString)
  RNAChars
end

function alphabet(protein::ProteinString)
  ProteinChars
end

# Frequency counts for the bases in a genetic sequence.
function countChars(gen::GeneticString)

  alpha = alphabet(gen)

  indices = [ alpha[i]=>i for i=1:length(alpha) ]

  function loop(counts::Array{Int, 1}, next::Char)
    counts[indices[next]] += 1
    counts
  end

  foldl(loop, zeros(Int, length(alpha)), gen.seq)
end

# Related problem - http://rosalind.info/problems/rna
function toRNAString(dna::DNAString)
  RNAString(replace(dna.seq, "T", "U"))
end

function toDNAString(rna::RNAString)
  DNAString(replace(rna.seq, "U", "T"))
end

function toProteinString(rna::RNAString)

  res = String[]

  for i = 1:3:length(rna)-2
    push!(res,RNACodonTable[rna.seq[i:i+2]])
  end
  ProteinString(join(res))
end

# Related problem - http://rosalind.info/problems/revc
function complement(dna::DNAString)

  complements = ['A'=>'T',
                 'C'=>'G',
                 'G'=>'C',
                 'T'=>'A']

  DNAString(map((ch) -> complements[ch], dna.seq))
end

function complement(rna::RNAString)

  complements = ['A'=>'U',
                 'C'=>'G',
                 'G'=>'C',
                 'U'=>'A']

  RNAString(map((ch) -> complements[ch], rna.seq))
end

function length(gen::GeneticString)
  length(gen.seq)
end

# Related problem - http://rosalind.info/problems/revc
function reverse(dna::DNAString)
  DNAString(reverse(dna.seq))
end

function reverse(rna::RNAString)
  RNAString(reverse(rna.seq))
end

function gcContent(dna::DNAString)

  GC = ['G', 'C']

  function loop(count::Int64, base::Char)
    if base in GC
      count + 1
    else
      count
    end
  end

  foldl(loop, 0, dna.seq) / length(dna)
end

# Related problem - http://rosalind.info/problems/hamm
function hammingDistance(seq1::GeneticString, seq2::GeneticString)

  seq1Len = length(seq1)
  seq2Len = length(seq2)

  # ensure sequence lengths are equal
  seq1Len == seq2Len ||
  error("Sequences must have same length.")

  function loop(count::Int64, index::Int64)
    if seq1.seq[index] != seq2.seq[index]
      count + 1
    else
      count
    end
  end

  # count mismatches of corresponding bases
  foldl(loop, 0, 1:seq1Len)
end

function findMotif(seq::String, motif::String)

  seqLen = length(seq)
  motifLen = length(motif)

  seqLen >= motifLen || error("motif length excedes sequence length")

  foldl((acc, i) ->
    if seq[i:i+motifLen-1] == motif
      push!(acc, i)
    else
      acc
    end, Int[], 1:seqLen - motifLen + 1)

end

function loadText(filename::String)

  # read text from file
  raw = open(readall, filename)

  replace(strip(raw), r"\n|\r", "")
end

function loadDNAString(filename::String)
  DNAString(loadText(filename))
end

function loadRNAString(filename::String)
  RNAString(loadText(filename))
end

function linesFromFile(filename::String)

  raw = open(readlines, filename)
  map(strip, raw)
end

end
