#' FastaUtils: Utilities for DNA/RNA sequence processing
#'
#' The package \pkg{FastaUtils} provides tools for the manipulation of DNA/RNA sequence data.
#'
#'@details The package \pkg{FastaUtils} depends on \pkg{Biostrings} which can be installed from Bioconductor at http://bioconductor.org/packages/release/bioc/html/Biostrings.html
#' 
#' To see the preferable citation of the package, type citation("FastaUtils").
#'@docType package
#'@name FastaUtils
#'@author Guillem Salazar <salazar@@icm.csic.es>

NULL


#' Select a subset of sequences from a fasta file
#'
#' This function loads a fasta file, selects a subset on sequences based on the sequence's names and writes a new fsata file.
#' @param file Fasta file
#' @param subset Vector with (exact) names of the sequences to be retrieved.
#' @param out Path to output file. If absent, the '.mod' is added to the input file's name (path is thus conserved).
#' @keywords FastaUtils
#' @return Writes a fasta file with the selected sequences.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' subset.fasta(file="http://greengenes.lbl.gov/Data/JD_Tutorial/UnAligSeq24606.txt",subset=c("24.6jsd1.Tut","24.6jsd2.Tut ","24.6jsd3.Tut "),out="out.fasta")

subset.fasta<-function(file=NULL,subset=NULL,out=paste(file,".subset",sep="")){
  library(Biostrings)
  sequences<-readDNAStringSet(file)
  if (all(as.character(subset) %in% names(sequences))==FALSE) stop("There are names in 'subset' not present in the fasta file")
  pos<-match(as.character(subset),names(sequences))
  writeXStringSet(sequences[pos],filepath=out)}

#' Convert a DNA fasta to an RNA fasta
#'
#' This function loads a fasta filecontaining DNA sequences and writes a fasta file with the same sequences converted to its RNA equivalent.
#' @param file Fasta file
#' @param out Path to output file. If absent, the '.rna' is added to the input file's name (path is thus conserved).
#' @keywords FastaUtils
#' @return Writes a fasta file.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' DNA2RNA(file="http://greengenes.lbl.gov/Data/JD_Tutorial/UnAligSeq24606.txt",out="out.fasta")

DNA2RNA<-function(file=NULL,out=paste(file,".rna",sep="")){
  library(Biostrings)
  sequences<-readDNAStringSet(file)
  writeXStringSet(RNAStringSet(sequences),filepath=out)}

#' Convert a RNA fasta to an DNA fasta
#'
#' This function loads a fasta filecontaining RNA sequences and writes a fasta file with the same sequences converted to its DNA equivalent.
#' @param file Fasta file
#' @param out Path to output file. If absent, the '.dna' is added to the input file's name (path is thus conserved).
#' @keywords FastaUtils
#' @return Writes a fasta file.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' # See 'DNA2RNA()' function

RNA2DNA<-function(file=NULL,out=paste(file,".dna",sep="")){
  library(Biostrings)
  sequences<-readRNAStringSet(file)
  writeXStringSet(DNAStringSet(sequences),filepath=out)}
