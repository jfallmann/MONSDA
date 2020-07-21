=================================================
The ICS (IdentifierConditionSetting relationship)
=================================================

For each ID you work on, you can define one or multiple conditions and
settings that will be used for the analysis). The ICS also sets the
file structure to follow for the FASTQ directory, where the ID is the
first level and the Condition the second. Setting is used by
``RunSnakemake.py|RunNextflow.py`` to enable processing of the same
samples under different settings like mapping tools, trimming tools
and later also postprocessing tools or commandline options for these
tools.

As an example, I want to analyse samples retreived from LabA on
01012020 (yes that happens), with the mapping tools star and segemehl,
my ICS would look like this
``LabA:01012020:star,LabA:01012020:segemehl`` and my FASTQ directory
would resemble that like ``FASTQ/LabA/01012020``. The '01012020'
directory would thereby contain all the fastq.gz files I need for
analysis. This works of course also if you want to analyze samples
from different dates and same lab with same settings or different labs
and so on.
