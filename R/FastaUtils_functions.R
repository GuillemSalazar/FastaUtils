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

#' Select a random sample of sequences
#'
#' This function select a random sample of sequences with or without replacement.
#' @param infile Multi-Fasta file
#' @param nseq Sample size.
#' @param file.out Path to output file.
#' @param replacement Logical for using or not (default) sampling with replacement.
#' @keywords FastaUtils
#' @return Writes a fasta file with sampled sequences.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' fasta.sample(infile="http://greengenes.lbl.gov/Data/JD_Tutorial/UnAligSeq24606.txt",nseq=10,file.out="out.fasta")

fasta.sample<-function(infile=NULL,nseq=NULL,file.out=NULL,replacement=FALSE){
  library(Biostrings)
  seqs<-readDNAStringSet(infile)
  selected<-seqs[sample(1:length(seqs),nseq,replace=replacement)]
  writeXStringSet(selected,filepath=file.out)}

#' Cuts peaces of a fixed length and a random position for all sequences or for a random sample of sequences
#'
#' This function cuts peaces of sequences of a fixed length ('read.length') either by sampling randomly a multi-fasta file or by selecting an equal number of peaces from each sequence within a multi-fasta file.
#' @param infile Multi-Fasta file.
#' @param sampling.type Sampling type: select a fixed peaces from each sequence ('none', default) or performs a uniform random sampling ('uniform'). In the first case the total number of peaces has to be multiple of the number of sequences.
#' @param total.reads Total number of peaces.
#' @param read.length Length of te peaces to be cut.
#' @param file.out Path to output file.
#' @param replacement Logical for using or not (default) sampling of sequences with replacement. Only applies if 'sampling.type=T'.
#' @keywords FastaUtils
#' @return Writes a fasta file with the selected sequences and a tab-delimites file with the sequence name and the start and end position of the cut. This file gets the 'file.out' names and adds '.info' at the end.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' fasta.cutter(infile="http://greengenes.lbl.gov/Data/JD_Tutorial/UnAligSeq24606.txt",sampling.type="uniform",file.out="out.fasta",total.reads=15,read.length=20)


fasta.cutter<-function(infile=NULL,sampling.type="none",total.reads=NULL,read.length=NULL,file.out=NULL,replacement=FALSE){
  library(Biostrings)
  seqs<-readDNAStringSet(infile)
  sampling.type<-match.arg(sampling.type,c("none","uniform"))
  reads.per.seq<-total.reads/length(seqs)
  
  if (sampling.type=="none" & reads.per.seq==as.integer(reads.per.seq)) positions<-rep(1:length(seqs),each=reads.per.seq)
  if (sampling.type=="none" & reads.per.seq!=as.integer(reads.per.seq)) stop("if sampling.type='none' total.reads should be a multiple of infile's number of sequences.","\ntotal.reads: ",total.reads,"\nnumber of sequences: ",length(seqs),"\nreads/sequences= ",reads.per.seq)
  if (sampling.type=="uniform") positions<-sample(1:length(seqs),total.reads,replace=replacement)
  
  
  seqs.new<-seqs[positions]
  seqs.new.length<-unlist(lapply(seqs.new,length))
  seqs.new.ini.max<-seqs.new.length-read.length
  get.ini<-function(x=NULL){sample(1:x,1)}
  seqs.new.ini<-unlist(lapply(seqs.new.ini.max,get.ini))
  seqs.new.end<-seqs.new.ini+read.length-1
  
  seqs.new<-subseq(seqs.new,start=seqs.new.ini,end=seqs.new.end)
  names(seqs.new)<-paste(names(seqs.new),seqs.new.ini,seqs.new.end,sep="_")
  writeXStringSet(seqs.new,filepath=file.out)
  write.table(data.frame(sequence=names(seqs.new),start.pos=seqs.new.ini,end.pos=seqs.new.end),sep="\t",col.names=NA,quote=F,file=paste(file.out,"info",sep="."))}
