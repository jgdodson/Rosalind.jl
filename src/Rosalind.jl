module Rosalind

import Base.length
import Base.reverse

# (✓) Create abstract type for genetic sequence
# (✓) Add concrete protein sequence type
# ( ) Finish candidates function
# ( ) Remove inefficient repeated concatenation from candidates function
# ( ) Finish working on toProteinString
# ( ) Compare speeds of direct iteration and foldl

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

const DNACodonTable = ["TTT" => "F", "CTT" => "L", "ATT" => "I",
                       "GTT" => "V", "TTC" => "F", "CTC" => "L",
                       "ATC" => "I", "GTC" => "V", "TTA" => "L",
                       "CTA" => "L", "ATA" => "I", "GTA" => "V",
                       "TTG" => "L", "CTG" => "L", "ATG" => "M",
                       "GTG" => "V", "TCT" => "S", "CCT" => "P",
                       "ACT" => "T", "GCT" => "A", "TCC" => "S",
                       "CCC" => "P", "ACC" => "T", "GCC" => "A",
                       "TCA" => "S", "CCA" => "P", "ACA" => "T",
                       "GCA" => "A", "TCG" => "S", "CCG" => "P",
                       "ACG" => "T", "GCG" => "A", "TAT" => "Y",
                       "CAT" => "H", "AAT" => "N", "GAT" => "D",
                       "TAC" => "Y", "CAC" => "H", "AAC" => "N",
                       "GAC" => "D", "TAA" => "X", "CAA" => "Q",
                       "AAA" => "K", "GAA" => "E", "TAG" => "X",
                       "CAG" => "Q", "AAG" => "K", "GAG" => "E",
                       "TGT" => "C", "CGT" => "R", "AGT" => "S",
                       "GGT" => "G", "TGC" => "C", "CGC" => "R",
                       "AGC" => "S", "GGC" => "G", "TGA" => "X",
                       "CGA" => "R", "AGA" => "R", "GGA" => "G",
                       "TGG" => "W", "CGG" => "R", "AGG" => "R",
                       "GGG" => "G"]

# Superclass of the concrete types in this module.
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

  const alpha = alphabet(gen)

  const indices = [ alpha[i]=>i for i=1:length(alpha) ]

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

# Related problem - http://rosalind.info/problems/prot
function toProteinString(rna::RNAString)

  res = String[]

  for i = 1:3:length(rna)-2
    push!(res,RNACodonTable[rna.seq[i:i+2]])
  end
  ProteinString(join(res))
end

function toProteinString(dna::DNAString)

  res = String[]

  for i = 1:3:length(dna)-2
    push!(res,DNACodonTable[dna.seq[i:i+2]])
  end
  ProteinString(join(res))
end


# Related problem - http://rosalind.info/problems/revc
function complement(dna::DNAString)

  const complements = ['A'=>'T',
                 'C'=>'G',
                 'G'=>'C',
                 'T'=>'A']

  DNAString(map((ch) -> complements[ch], dna.seq))
end

function complement(rna::RNAString)

  const complements = ['A'=>'U',
                 'C'=>'G',
                 'G'=>'C',
                 'U'=>'A']

  RNAString(map((ch) -> complements[ch], rna.seq))
end

function length(gen::GeneticString)
  length(gen.seq)
end

function reverse(dna::DNAString)
  DNAString(reverse(dna.seq))
end

function reverse(rna::RNAString)
  RNAString(reverse(rna.seq))
end

# Related problem - http://rosalind.info/problems/revc
function reverseComplement(dnaOrRna::Union(DNAString, RNAString))
  reverse(complement(dnaOrRna))
end

# Related problem - http://rosalind.info/problems/gc
function gcContent(dnaOrRna::Union(DNAString, RNAString))

  const GC = ['G', 'C']

  function loop(count::Int64, base::Char)
    if base in GC
      count + 1
    else
      count
    end
  end

  foldl(loop, 0, dnaOrRna.seq) / length(dnaOrRna)
end

# Related problem - http://rosalind.info/problems/hamm
function hammingDistance(seq1::GeneticString, seq2::GeneticString)

  const seq1Len = length(seq1)
  const seq2Len = length(seq2)

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

function findMotif(gen::GeneticString, motif::String)

  const seqLen = length(seq)
  const motifLen = length(motif)

  seqLen >= motifLen || error("motif length excedes sequence length")

  function loop(acc::Array{Int64, 1}, index::Int64)
    if gen.seq[index:index+motifLen-1] == motif
      push!(acc, index)
    else
      acc
    end
  end

  foldl(loop, Int[], 1:seqLen - motifLen + 1)

end

# Related problem - http://rosalind.info/problems/orf
function candidates(dna::DNAString)

  function helper(dna::DNAString)
    startCodon = "ATG"
    stopCodons = ["TAG", "TGA", "TAA"]

    bins = [i=>String[] for i=0:2]
    acc = String[]

    for i=1:length(dna)-2
      codon = dna.seq[i:i+2]
      binNum = i % 3
      if codon == startCodon
        map!((b) -> string(b, DNACodonTable[codon]), bins[binNum])
        push!(bins[binNum], DNACodonTable[startCodon])
      elseif codon in stopCodons
        append!(acc, bins[binNum])
        bins[binNum] = String[]
      else
        map!((b) -> string(b, DNACodonTable[codon]), bins[binNum])
        end
    end

    acc
  end

  append!(helper(dna), helper(reverseComplement(dna)))
end

function candidates2(dna::DNAString)

  startCodon = "ATG"
  stopCodons = ["TAG", "TGA", "TAA"]

  bins = [Array{String, 1}[] for i=1:3]
  acc = Array{String, 1}[]

  for i=1:length(dna)-2
    codon = dna.seq[i:i+2]
    binNum = 1 + (i % 3)
    if codon == startCodon
      map!((b) -> push!(b, DNACodonTable[codon]), bins[binNum])
      push!(bins[binNum], [DNACodonTable[startCodon]])
    elseif codon in stopCodons
      append!(acc, bins[binNum])
      bins[binNum] = String[]
    else
      map!((b) -> push!(b, DNACodonTable[codon]), bins[binNum])
    end
  end

  map(join, acc)
end

function candidates3(dna::DNAString)

  startCodon = "ATG"
  stopCodons = ["TAG", "TGA", "TAA"]

  bins = [Array{Int, 1}[] for i=1:3]
  acc = Array{Int, 1}[]

  for i=1:length(dna)-2
    codon = dna.seq[i:i+2]
    binNum = 1 + (i % 3)
    if codon == startCodon
      for item in bins[binNum]
        item[2] += 1
      end
      push!(bins[binNum], [i, i])
    elseif codon in stopCodons
      append!(acc, bins[binNum])
      bins[binNum] = Array{Int, 1}[]
    else
      for item in bins[binNum]
        item[2] += 1
      end
    end
  end

  acc
end

function candidates4(dna::DNAString)

  startCodon = "ATG"
  stopCodons = ["TAG", "TGA", "TAA"]

  bins = [ i=> Array(Int, 0) for i=0:2]
  acc = Array(Int, 0)

  for i=1:length(dna)-2
    codon = dna.seq[i:i+2]
    binNum = i % 3
    if codon == startCodon
      for j in 2:2:length(bins[binNum])
        bins[binNum][j] += 1
      end
      push!(bins[binNum], i, i)
    elseif codon in stopCodons
      append!(acc, bins[binNum])
      bins[binNum] = Array(Int, 0)
    else
      for j in 2:2:length(bins[binNum])
        bins[binNum][j] += 1
      end
    end
  end

  acc
end

function loadText(filename::String)

  open(readall, filename)
end

function cleanText(raw::String)

  replace(strip(raw), r"\n|\r", "")
end

function readFastaFile(filename::String)

  fastaPattern = r">([^\s]*) ?([^\n]*)\n([^>]*)"s
  raw = loadText(filename)

  map((m) -> DNAString(cleanText(m.captures[3])), eachmatch(fastaPattern, raw))
end

function loadDNAString(filename::String)
  DNAString(cleanText(loadText(filename)))
end

function loadRNAString(filename::String)
  RNAString(cleanText(loadText(filename)))
end

function loadProteinString(filename::String)
  ProteinString(cleanText(loadText(filename)))
end

function linesFromFile(filename::String)

  raw = open(readlines, filename)
  map(strip, raw)
end

end

