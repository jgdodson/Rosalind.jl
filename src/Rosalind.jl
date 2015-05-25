module Rosalind

# (✓) Create abstract type for genetic sequence
# ( ) Add concrete protein sequence type

abstract GenSeq


immutable DNASeq <: GenSeq
  seq::String

  DNASeq(seq) = isDNAString(seq) ? new(seq) : error("String contains non-ACTG character.")
end

immutable RNASeq <: GenSeq
  seq::String

  RNASeq(seq) = isRNAString(seq) ? new(seq) : error("String contains non-ACGU character.")
end

function isDNAString(seq::String)

  bases = ['A','C','G','T']

  all(b -> b in bases, seq)
end

function isRNAString(seq::String)

  bases = ['A','C','G','U']

  all(b -> b in bases, seq)
end


# Frequency counts for the bases in a DNA sequence.
function countBasesDNA(dna::DNASeq)

  counts::Array{Int64, 1} = zeros(4)

  for c in dna.seq
    if c == 'A'
      counts[1] += 1
    elseif c == 'C'
      counts[2] += 1
    elseif c == 'G'
      counts[3] += 1
    else
      counts[4] += 1
    end
  end

  counts
end

# Related problem - http://rosalind.info/problems/rna
function transcribe(dna::DNASeq)
  RNASeq(replace(dna.seq, "T", "U"))
end

# Related problem - http://rosalind.info/problems/revc
function complement(dna::DNASeq)

  complements = ['A'=>'T',
                 'C'=>'G',
                 'G'=>'C',
                 'T'=>'A']

  DNASeq(map((ch) -> complements[ch], dna.seq))
end

function complement(rna::RNASeq)

  complements = ['A'=>'U',
                 'C'=>'G',
                 'G'=>'C',
                 'U'=>'A']

  RNASeq(map((ch) -> complements[ch], rna.seq))
end

function length(seq::GenSeq)
  Base.length(seq.seq)
end

# Related problem - http://rosalind.info/problems/revc
function reverse(dna::DNASeq)
  DNASeq(reverse(dna.seq))
end

function reverse(rna::RNASeq)
  RNASeq(reverse(rna.seq))
end

# Related problem - http://rosalind.info/problems/hamm
function hammingDistance(seq1::GenSeq, seq2::GenSeq)

  seq1Len = length(seq1)
  seq2Len = length(seq2)

  # ensure sequence lengths are equal
  seq1Len == seq2Len ||
  error("Sequences must have same length.")

  # count mismatches of corresponding bases
  foldl((acc, i) ->
    if seq1.seq[i] != seq2.seq[i]
      acc + 1
    else
      acc
    end, 0, 1:seq1Len)
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

function loadDNASeq(filename :: String)

  # read text from file
  raw = open(readall, filename)

  clean = replace(strip(raw), r"\n|\r", "")

  DNASeq(clean)
end

function linesFromFile(filename::String)

  raw = open(readlines, filename)
  map(strip, raw)
end

end
