from .demultiplex_utils import demultiplex
# from .bowtie2_utils import bowtie2
from .bam2dinucleotide_utils import bam2dinucleotide
from .statistic_utils import statistic
from .report_utils import report
from .fig_genome_utils import fig_genome
from .fig_position_utils import fig_position
from .common_utils import check_file_exists, get_reference_start, runCmd, readFileMap, checkFileMap, remove_chr, write_r_script

from .__version__ import __version__
