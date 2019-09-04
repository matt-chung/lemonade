#!/usr/bin/env sh

VER=1.0.0
PROJECT_NAME="ggplot_sashimi.sh"

R_BIN_DIR="/home/mattchung/.local/bin/R-3.4.4/bin"
BAMTOOLS_BIN_DIR="/usr/local/packages/bamtools-2.5.1/bin"
REGTOOLS_BIN_DIR="/usr/local/packages/regtools-0.5.0/bin"
SAMTOOLS_BIN_DIR="/usr/local/packages/samtools-1.9/bin"


showhelp() {
#  _initpath
  version
  echo "Usage: ggplot_sashimi.sh ...[parameters]....
  Commands:
  --help, -h               Show this help message.
  --version, -v            Show version info.
  --bam, -i                Path to BAM file
  --coords , -c            Specifies location to generate a sashimi plot (CONTIG:STARTPOS-ENDPOS)
  --stranded, -s           Indicates if reads are from a strand-specific assay <yes, no, reverse> (default: no)
  --output-dir, -o         Path to output directory
  --threads, -@            Number of threads to use for BAM sorting and indexing
  "
}

sortbam(){
  if [ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" ]; then
      "$SAMTOOLS_BIN_DIR"/samtools sort "$BAM" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")"
  fi
}

splitbam(){
  if [ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.f.bam/")" && "$STRANDED" = "yes" ]; then
		"$BAMTOOLS_BIN_DIR"/bamtools filter -in "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" -isReverseStrand FALSE -isMateReverseStrand FALSE -out "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.f.bam/")"
		"$BAMTOOLS_BIN_DIR"/bamtools filter -in "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" -isReverseStrand TRUE -isMateReverseStrand TRUE -out "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.r.bam/")"
  elif [ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.f.bam/")" && "$STRANDED" = "reverse" ]; then
  	"$BAMTOOLS_BIN_DIR"/bamtools filter -in "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" -isReverseStrand TRUE -isMateReverseStrand TRUE -out "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.f.bam/")"
  	"$BAMTOOLS_BIN_DIR"/bamtools filter -in "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" -isReverseStrand FALSE -isMateReverseStrand FALSE -out "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.r.bam/")"
  fi
}

indexbam(){
  if [ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam.bai/")" & "$STRANDED" = "no" ]; then
  	"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")"
  elif [ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.f.bam.bai/")" & "$STRANDED" != "no" ]; then
  	"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.f.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.r.bam/")"
  fi
}

getplotcoords(){
	CONTIG="$(echo "$COORDS)" | cut -d":" -f1)"
	START="$(echo "$COORDS)" | cut -d":" -f1 | cut -d"-" -f1)"
	STOP="$(echo "$COORDS)" | cut -d":" -f1 | cut -d"-" -f2)"
	
	PLOT_START=$(($variableA/1000*1024)) 
	PLOT_STOP



}








of(







}

findjunctions(){
"$REGTOOLS_BIN_DIR"/regtools junctions extract "$BAM" -r -o "$BAM".junctions.tsv
}



main(){



}




