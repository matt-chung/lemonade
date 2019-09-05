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

splitbam(){
  if [[ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.f.bam/")" && "$STRANDED" != "no" ]]; then
		"$SAMTOOLS_BIN_DIR"/samtools view -bh -f 83 "$BAM" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.83.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools view -bh -f 99 "$BAM" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.99.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools view -bh -f 147 "$BAM" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.147.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools view -bh -f 163 "$BAM" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.163.bam/")"
		if [[ "$STRANDED" = "yes" ]]; then
			"$SAMTOOLS_BIN_DIR"/samtools merge -f -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.f.bam/")" "$( echo "$BAM" | sed "s/[.]bam$/.99.bam/g")" "$( echo "$BAM" | sed "s/[.]bam$/.147.bam/g")"
			"$SAMTOOLS_BIN_DIR"/samtools merge -f -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.r.bam/")" "$( echo "$BAM" | sed "s/[.]bam$/.83.bam/g")" "$( echo "$BAM" | sed "s/[.]bam$/.163.bam/g")"
		elif [[ "$STRANDED" = "reverse" ]]; then
			"$SAMTOOLS_BIN_DIR"/samtools merge -f -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.f.bam/")" "$( echo "$BAM" | sed "s/[.]bam$/.83.bam/g")" "$( echo "$BAM" | sed "s/[.]bam$/.163.bam/g")"
			"$SAMTOOLS_BIN_DIR"/samtools merge -f -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.r.bam/")" "$( echo "$BAM" | sed "s/[.]bam$/.99.bam/g")" "$( echo "$BAM" | sed "s/[.]bam$/.147.bam/g")"
		fi
		rm "$( echo "$BAM" | sed "s/[.]bam$/.99.bam/g")" 
		rm "$( echo "$BAM" | sed "s/[.]bam$/.147.bam/g")"
		rm "$( echo "$BAM" | sed "s/[.]bam$/.83.bam/g")"
		rm "$( echo "$BAM" | sed "s/[.]bam$/.163.bam/g")"
	fi
}

sortbam(){
	if [[ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.f.sortedbyposition.bam/")" && "$(echo "$BAM" | sed "s/[.]bam$/.r.sortedbyposition.bam/")" && "$STRANDED" != "no" ]]; then
		"$SAMTOOLS_BIN_DIR"/samtools sort "$(echo "$BAM" | sed "s/[.]bam$/.f.bam/")" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.f.sortedbyposition.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools sort "$(echo "$BAM" | sed "s/[.]bam$/.r.bam/")" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.r.sortedbyposition.bam/")"
	elif [[ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" ]]
		"$SAMTOOLS_BIN_DIR"/samtools sort "$BAM" -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")"
	fi
}

indexbam(){
	if [[ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.f.sortedbyposition.bam.bai/")" && ! -f "$(echo "$BAM" | sed "s/[.]bam$/.r.sortedbyposition.bam.bai/")" && "$STRANDED" != "no" ]]; then
		"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.f.sortedbyposition.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.r.sortedbyposition.bam/")"
	elif [[ ! -f "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam.bai/")" ]]; then
		"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")"
	fi
}

getplotcoords(){
	CONTIG="$(echo "$COORDS" | cut -d":" -f1)"
	START="$(echo "$COORDS" | cut -d":" -f2 | cut -d"-" -f1)"
	STOP="$(echo "$COORDS" | cut -d":" -f2 | cut -d"-" -f2)"

	PLOT_START=$(($START - ($STOP - $START + 1)/10))
	PLOT_STOP=$(($STOP + ($STOP - $START + 1)/10))
}

getjunctions(){
	if [ "$STRANDED" = "no" ]; then
		"$REGTOOLS_BIN_DIR"/regtools junctions extract "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" -r "$CONTIG":"$PLOT_START"-"$PLOT_STOP" -o "$OUTPUT-DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".junctions.tsv/")"
	else
		"$REGTOOLS_BIN_DIR"/regtools junctions extract "$(echo "$BAM" | sed "s/[.]bam$/.f.sortedbyposition.bam/")" -r "$CONTIG":"$PLOT_START"-"$PLOT_STOP" -o "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".f.junctions.tsv/")"
		"$REGTOOLS_BIN_DIR"/regtools junctions extract "$(echo "$BAM" | sed "s/[.]bam$/.r.sortedbyposition.bam/")" -r "$CONTIG":"$PLOT_START"-"$PLOT_STOP" -o "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".r.junctions.tsv/")"
	fi
}

getmpileup(){
	if [ "$STRANDED" = "no" ]; then
		"$SAMTOOLS_BIN_DIR"/samtools mpileup -a -d 1000000 -r "$CONTIG":"$PLOT_START"-"$PLOT_STOP" -o "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/.mpileup/")" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/")" 
	else
		"$SAMTOOLS_BIN_DIR"/samtools mpileup -a -d 1000000 -r "$CONTIG":"$PLOT_START"-"$PLOT_STOP" -o "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".f.mpileup/")" "$(echo "$BAM" | sed "s/[.]bam$/.f.sortedbyposition.bam/")"
		"$SAMTOOLS_BIN_DIR"/samtools mpileup -a -d 1000000 -r "$CONTIG":"$PLOT_START"-"$PLOT_STOP" -o "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".r.mpileup/")" "$(echo "$BAM" | sed "s/[.]bam$/.r.sortedbyposition.bam/")"
	fi
}

summarizempileup(){
	echo -e ""contig"\t"position"\t"depth"\t"A"\t"T"\t"G"\t"C"\t"a"\t"t"\t"g"\t"c"\t"insertions"\t"deletions"" > "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".f.mpileup.summary/")"
	
	while read LINE; do
		CONTIG="$(echo "$LINE" | cut -f1)"
		POSITION="$(echo "$LINE" | cut -f2)"
		REF_BASE="$(echo "$LINE" | cut -f3)"
		REV_REF_BASE="$(echo "$REF_BASE" | tr '[:upper:]' '[:lower:]')"
		BASES="$(echo "$LINE" | cut -f4 | sed "s/\^.//g" | sed "s/[<>$]//g" | sed "s/[.]/"$REF_BASE"/g" | sed "s/[,]/"$REV_REF_BASE"/g")"
		INSERTION=""
		DELETION=""
		
		BASES="$(echo "$BASES" | sed "s/[<>$]//g")"
		while [ "$(echo $BASES | grep -ob "[+-]" | head -n1 | wc -l)" -gt 0 ]; do
			INDEL_START="$(echo $BASES | grep -ob "[+-]" | head -n1 | cut -d":" -f1)"
			INDEL_LENGTH="${BASES:${INDEL_START} + 1:${INDEL_START} + 1}"
			INDEL="${BASES:${INDEL_START}:${INS_LENGTH} + 2}"
			OCCURENCES="$(echo "$BASES" | grep -o -- "$INDEL" | wc -l)"
			if [ "$(echo $INDEL | grep -ob "+" | wc -l)" -gt 0 ]; then
				INSERTION=""$INSERTION"|"$OCCURENCES":"$INDEL""
			else
				DELETION=""$DELETION"|"$OCCURENCES":"$INDEL""
			fi
			BASES="$(echo "$BASES" | sed "s/"$INDEL"//g")"
		done

		A="$(echo "$BASES" | grep -o "A" | wc -l)
		T="$(echo "$BASES" | grep -o "T" | wc -l)
		G="$(echo "$BASES" | grep -o "G" | wc -l)
		C="$(echo "$BASES" | grep -o "C" | wc -l)
		
		a="$(echo "$BASES" | grep -o "a" | wc -l)
		t="$(echo "$BASES" | grep -o "t" | wc -l)
		g="$(echo "$BASES" | grep -o "g" | wc -l)
		c="$(echo "$BASES" | grep -o "c" | wc -l)
		
		DEPTH="$(echo "$LINE" | wc -c)"
		
		echo -e ""$CONTIG"\t"$POSITION"\t"$DEPTH"\t"$A"\t"$T"\t"$G"\t"$C"\t"$a"\t"$t"\t"$g"\t"$c"\t"$INSERTION"\T"$DELETION"" >>  "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".f.mpileup.summary/")"
	done <  "$OUTPUT_DIR"/"$(basename "$BAM" | sed "s/[.]bam$/."$CONTIG":"$PLOT_START"-"$PLOT_STOP".f.mpileup/")"
}

main(){
	splitbam
	sortbam
	indexbam
	getjunctions
	getmpileup
	summarizempileup
	
	createRinputs
}

main
