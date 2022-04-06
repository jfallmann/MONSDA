.. _preparation:

Plannig a project run
======================

Before you set up a ``MONSDA`` project, you should make a few preliminary considerations. For example, if you want to process data from an experiment, you need to make sure that ``MONSDA`` "understands" your experiment. 
You should also have a good overview of your data and know in advance how exactly it should be calculated.

To be well prepared, go through the following steps. Afterwards you will run the configurator, which helps you to set up ``MONSDA``.

STEP 1: What do i want to do exactly?
-------------------------------------

  * What is the experimental design? Which conditions do i have?

    If you have data from an experiment and it is necessary to differ bewtween conditions, ``MONSDA`` represents it as a ``condition-tree``. 
    Make sure, you know how to transfer your Experimental design to the ``condition-tree``. An example could be a simple 'knockout' Experiment. 
    Your would transfer it as following: 
    
    .. code-block::

                            ┌──Condition1: "Wildtype" 
      Knockout-Experiment ──┤
                            └──Condition2: "Knockout"

    Your experimental design can be more komplex. To understand better the concept of the ``condition-tree`` click :ref:`here<condition-tree>`.

  * What kind of analysis i want to do? Which processing steps and tools i need?

    You probably have an idea of the workflow you want to create. You should write down which processing steps you will establish. 
    ``MONSDA`` already offers a lot of processing tools. Jump to the caption :ref:`Workflow Overview<WFoverview>` to get out 
    which tools you want to use for your workflow.

  * Which settings should the tools run with?

    Depending on your data and your specific analysis most of the tools must be assigned special settings or options. 
    Of course it is possible to adjust the settings afterwards, but it is usefull to think about it beforehand wich settings you want to run each tool with.
    For example if you integrate a mapper tool to your Workflow, you should know which option-flags an settings you need. 
    It can be helpfull to test your settings by hand before you integrate them in the ``MONSDA`` configuration. 
    To understand better how ``MONSDA`` integrates tool options have a look at :ref:`the config-file<config-file>`.

STEP 2: Keep the overview
-------------------------

  You may want to do a common sequencing analysis. Most likely you will need a genome fasta file and annotation files in addition to your samples.
  Maybe there are further additional files your Workflow need.
  
  So before you start the configurator and set up the ``MONSDA`` project, make sure you have organized your data well. The configurator will establish a project 
  folder with softlinks to your data. It would be bad if you move your data afterwards.   


STEP 3: Hardware issues
-----------------------

  ``MONSDA`` can be run on different platforms and also allows to controll multithreading. If you want to run it on Slurm, have a look :ref:`here<slurm>`. 
  Also make sure how many threads you can use on your machine.
  
  Last but not least the processing results will be saved in your project-folder. This can take up a lot of disk space. make sure you have enough left over. 
  Another question is on which pc you want to start the pipeline


