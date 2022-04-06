.. _preparation:

Planning a project run
======================

Before you set up a ``MONSDA`` project, you should take a few things into considerations to make sure that ``MONSDA`` "understands" your experiment. 
You should also have a good overview of your data and best decide in advance what analysis steps you want to run. Although it is possible to run ``MONSDA`` sequentially, building up a config step by step, this will lead to some workflows being run multiple times, as a unique "tool-key" will be defined for the set of enabled workflow steps and configured tools, more information on that can be found in :ref:`here<condition-tree>`.

To be well prepared, go through the following steps. Afterwards you will run the ``Configurator``, which helps you to set up ``MONSDA``.

STEP 1: What do I want to do exactly?
-------------------------------------

  * What is the experimental design? Which conditions do I have?

    If you have data from an experiment and it is necessary to differ between conditions, ``MONSDA`` represents it as a ``condition-tree``. 
    Make sure, you know how to transfer your Experimental design to the ``condition-tree``. An example could be a simple 'knockout' Experiment. 
    Your would transfer it as following: 
    
    .. code-block::

                            ┌──Condition1: "Wild-Type" 
      Knockout-Experiment ──┤
                            └──Condition2: "Knockout"

    Your experimental design can be more complex. To get a better understanding of the concept behind the ``condition-tree`` click :ref:`here<condition-tree>`.

  * What kind of analysis do I want to run? Which processing steps and tools do I need?

    You probably have an idea of the workflow you want to create. You should write down which processing steps you will establish. 
    ``MONSDA`` already offers a bunch of processing tools. Jump to :ref:`Workflow Overview<WFoverview>` to get an overview of 
    which tools you can use for your workflow.

  * Which settings should the tools run with?

    Depending on your data and your specific analysis most of the tools must be assigned special settings or options. 
    Of course it is possible to adjust the settings afterwards, but this requires to delete generated results and is more tedious then to think about settings beforehand. ``MONSDA`` does not recommend settings for tools and will run all tools with their default settings. The only exception from that rule is a "paired" or "stranded" setting, which will automatically be derived from the config and resolved for each tool, where applicable.
    
    It can be helpful to test your settings by hand before you integrate them in the ``MONSDA`` configuration, however, you can always re-run ``MONSDA`` after cleaning up results and changing the configuration.
    To understand better how ``MONSDA`` integrates tool options have a look at :ref:`the config-file<config-file>`.

STEP 2: Keep the overview
-------------------------

  You may want to do a common sequencing analysis. Most likely you will need a genome FASTA file and annotation files in addition to your samples.
  Maybe there are further additional files your workflow needs.
  
  So before you start the configurator and set up the ``MONSDA`` project, make sure you have organized your data well. The configurator will establish a project 
  folder with softlinks to your data. Later changes require you to move your data afterwards, which is extra work that is not necessary.

If you later on share your workflow and configuration, we also assume that you share all input data including genome and annotation files, at least with links to the original download location, so it makes sense to keep that information.

STEP 3: Hardware issues
-----------------------

  ``MONSDA`` can be run on different platforms and allows to control multithreading. If you want to run it on *SLURM*, have a look :ref:`here<slurm>`. Always make sure you know how many threads you can run on your machine and keep an eye on memory consumption. For some workflow steps like mapping, memory can be a bottle-neck, so best check your *SLURM* profile and set caps on memory usage.
  
  Last but not least the processing results will be saved in your project-folder. This can take up a lot of disk space. make sure you have enough free space. 
  
  Details on input and output formats and files can be found in the :ref:`Workflow Overview<WFoverview>`.


