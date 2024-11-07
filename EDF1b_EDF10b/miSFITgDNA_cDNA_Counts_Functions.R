#A function for checking the frequency of a user specified list of sequences in a subseqset
Frequencies <- function(miSFITs, cgDNAs)
{
  #initialize a vector to store how frequently in the subseqset each variant in the validationset appears 
  Freqs = vector(length = length(miSFITs))
  totalreadcount = length(cgDNAs)
  for(i in 1:length(Freqs))
  {
    #reports progress to the user
    print(paste(Sys.time(), "progress is" ,i, "out of",length(miSFITs)))
    #grabs the next sequence in the user specified set
    sequence_to_count <- miSFITs[[i]]
    #counts that sequence
    Freqs[i] <- sum(vcountPattern(sequence_to_count,cgDNAs))/ totalreadcount
  }
  return(Freqs)
}
#A function for checking the absolute count of reads for of a user specified list of sequences in a subseqset
#This function differs from validation_variant_counter (which returns frequency), this function returns absolute read counts
Absolutes <- function(miSFITs, cgDNAs)
{
  #initialize a vector to store how frequently in the subseqset each variant in the validationset appears 
  Absols = vector(length = length(miSFITs))
  for(i in 1:length(Absols))
  {
    #reports progress to the user
    print(paste(Sys.time(), "progress is" ,i, "out of",length(Absols)))
    #grabs the next sequence in the user specified set
    sequence_to_count <- miSFITs[[i]]
    #counts that sequence
    Absols[i] <- sum(vcountPattern(sequence_to_count,cgDNAs))
  }
  return(Absols)
}
