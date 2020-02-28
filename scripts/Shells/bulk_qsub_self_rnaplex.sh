#RNAplex interactions cluster submitter
#!/bin/bash
#arguments $1 = path_fasta_directory, $2 = path_accessibility_directory $3 = path_output_directory
#Example call:
#qsub -t 1-2 ~egg/current/RNAnet/scripts/bulk_qsub_self_rnaplex.sh /scr/k70san/fall/RINTERwahn/PARIS/Annotated/Beds/SplitFasta/ /scr/k70san/fall/RINTERwahn/PARIS/Annotated/Beds/PlFold_Indexed/ /scr/k70san/fall/RINTERwahn/PARIS/PlexOut/

#$ -j yes
#$ -o /scr/k70san/fall/RINTERwahn/PARIS/PlexOut/RINTERwahn.out
#$ -e /scr/k70san/fall/RINTERwahn/PARIS/PlexOut/RINTERwahn.out
##$ -q hostname="xc00|xc01|xc02|xc03|xc05|xc06|xc07|xc08|tc00|tc01|tc03|tc02|tc04"
#$ -q short.q
#$ -N Muhhh
#BE CAREFUL ONLY ONE LINE FASTA
#iabulk_rnaplex args: path_fasta_directory fst_fasta_index snd_fasta_index path_accessibility_directory path_output_directory
/scr/k70san/fall/RINTERwahn/ia_bulk_rnaplex.sh $SGE_TASK_ID $SGE_TASK_ID $1 $2 $3
