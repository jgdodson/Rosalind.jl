module Rosalind

# Frequency counts for the bases in a DNA sequence.
function countBasesDNA(seq :: String)
  
  counts::Array{Int64, 1} = zeros(4)
  
  for c in seq
    if c == 'A'
      counts[1] += 1
    elseif c == 'C'
      counts[2] += 1
    elseif c == 'G'
      counts[3] += 1
    elseif c == 'T' 
      counts[4] += 1
    else error("DNA sequence contains non-ACGT character")
    end
  end
  
  counts
end

function transcribe(dna::String)
  replace(dna, "T", "U")
end

function complementDNA(dna :: String)

  map((ch) ->
    if ch == 'A' 
      'T'
    elseif ch == 'C' 
      'G'
    elseif ch == 'G' 
      'C'
    elseif ch == 'T' 
      'A'
    else
      error("DNA sequence contains non-ACTG character")
    end, dna)
end

function reverseComplementDNA(dna :: String)
  reverse(complementDNA(dna))
end

function countMutations(seq1::String, seq2::String)

  # ensure sequence lengths are equal
  length(seq1) == length(seq2) || 
  error("sequences must have same length")

  # count mismatches of corresponding bases
  foldl((acc, next) ->
    if seq1[next] != seq2[next]
      acc + 1
    else
      acc
    end, 0, 1:length(seq1))
end

function findMotif(seq :: String, motif :: String)

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

function loadSeq(filename :: String)
  
  # open datafile
  raw = open(readall, filename)
  
  replace(strip(raw), r"\n|\r", "")
end

function linesFromFile(filename::String)

  raw = open(readlines, filename)
  map(strip, raw)
end

end
