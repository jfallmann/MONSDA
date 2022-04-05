.. _condition-tree:

The Condition-Tree
==================

A key concept behind ``MONSDA`` is that config files are split according to conditions defined in the ``config`` file and each subworkflow is then run consecutively. This separates workflows in independent subworkflows, each potentially running on different input, with differing options and executables, without danger of interferring with each other.

As described in :ref:`preparation`, you should have an idea of what to analyse and how to split up your conditions if needed. For each ID you work on, you can define no, one or multiple conditions and settings that will be used for the analysis. The condition-tree also defines the directory structure to follow for Input and Output directories. We asume a general tree to look something like

.. code-block:: json

    'ID' -> 'CONDITION' -> 'SETTING'

Here *ID* is the first level and optional *Condition* the second. *Setting* is used by ``MONSDA`` to enable processing of the same samples under different settings or commandline options for e.g. mapping tools, trimming tools and later also postprocessing tools. ``MONSDA`` will also build an output directory based on the combination of tools used, e.g. fastqc-cutadapt-star-umitools, to indicate which combination of tools was used to generate the output and to prevent results from being mixed or overwritten.

As an example, I want to analyse samples retreived from LabA on 01012020 (yes that happens), with the mapping tools star and hisat, my condition-tree would look like this ``LabA:01012020`` and my FASTQ input directory would resemble that like ``FASTQ/LabA/01012020``. The '01012020' directory would thereby contain all the fastq.gz files I need for analysis as stated in the corresponding config file. As we assume that settings may change but the input files will stay the same, ``MONSDA`` will search one level upwards of the deepest level in the condition-tree. This means, if you have a tree like:

.. code-block:: json

    'ID1' -> 'CONDITION1' -> 'SETTING1', 'SETTING2', 'SETTING3'

You do not need to copy input from ``FASTQ/LabA/01012020`` to ``FASTQ/LabA/01012020/SETTING1/2/3``, instead ``MONSDA`` will find the input in ``FASTQ/LabA/01012020`` and generate output directories which contain the *Setting* level.

This works of course also if you want to analyze samples from different dates and same lab with same settings or different labs and so on. ```MONSDA``` will automatically generate output folders ``FASTQ/LabA/01012020/star`` and ``FASTQ/LabA/01012020/hisat`` if no other tools where configured to be used.

Optionally a user can also run one or the other tool with different settings for example to benchmark tools, e.g. ``map_stringent`` and ``map_relaxed`` and indicate this on the *Setting* level in the config file. FASTQ input will still be found in ``FASTQ/LabA/01012020`` while output files will appear in ``FASTQ/LabA/01012020/map_stringent/star`` and ``FASTQ/LabA/01012020/map_stringent/hisat`` or ``FASTQ/LabA/01012020/map_relaxed/star`` and ``FASTQ/LabA/01012020/map_relaxed/star`` respectively.
