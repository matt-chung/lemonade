#!/usr/bin/env sh

VER=1.0.0
PROJECT_NAME="ggplot_sashimi.sh"

R_BIN_DIR="/home/mattchung/.local/bin/R-3.4.4/bin"
REGTOOLS_BIN_DIR="/usr/local/packages/regtools-0.5.0/bin"
SAMTOOLS_BIN_DIR="/usr/local/packages/samtools-1.9/bin"


showhelp() {
#  _initpath
  version
  echo "Usage: ggplot_sashimi ...[parameters]....
  Commands:
  --help, -h               Show this help message.
  --version, -v            Show version info.
  --bam, -i                Path to BAM file
  --stranded, -s           Indicates if reads are from a strand-specific assay <yes, no, reverse> (default: no)
  "
}



