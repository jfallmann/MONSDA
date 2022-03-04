=================================================
The ICS (Identifier-Condition-Setting relationship)
=================================================

For each ID you work on, you can define no, one or multiple conditions and
settings that will be used for the analysis). The ICS also sets the
file structure to follow for the FASTQ input and all of the output directories, where the ID is the
first level and optional Condition the second. Setting is used by
```MONSDA``` to enable processing of the same
samples under different settings for e.g. mapping tools, trimming tools
and later also postprocessing tools or commandline options for these
tools. ```MONSDA``` will also build an output directory based on the combination of tools used,
e.g. fastqc-cutadapt-star-umitools, to indicate which combination of tools was used to
generate the output.

As an example, I want to analyse samples retreived from LabA on
01012020 (yes that happens), with the mapping tools star and hisat,
my ICS would look like this
``LabA:01012020`` and my FASTQ input directory
would resemble that like ``FASTQ/LabA/01012020``. The '01012020'
directory would thereby contain all the fastq.gz files I need for
analysis as stated in the corresponding config file. 
This works of course also if you want to analyze samples
from different dates and same lab with same settings or different labs
and so on. ```MONSDA``` will automatically generate output folders 
``FASTQ/LabA/01012020/star`` and ``FASTQ/LabA/01012020/hisat`` if no other tools
where configured to be used. 

Optionally a user can also run one or the other tool 
with different settings for example to benchmark tools,
e.g. ``map_stringent`` and ``map_relaxed`` and indicate this as Setting 
in the config file. FASTQ input will still be found in ``FASTQ/LabA/01012020``
while output files will appear in ``FASTQ/LabA/01012020/map_stringent/star`` and ``FASTQ/LabA/01012020/map_stringent/hisat`` or ``FASTQ/LabA/01012020/map_relaxed/star`` and ``FASTQ/LabA/01012020/map_relaxed/star`` respectively.