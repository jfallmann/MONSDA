.. role::  raw-html(raw)
    :format: html

================
Configure MONSDA
================


| MONSDA operates from a project folder containing a ``config`` json file.
| To enable easy setup and configuration, we host the executable ``Configurator``.
| The Configurator works as an interactive terminal user interface.
| It supports you in three main things:

1. Creating the complete initial project folder
2. Creating the configuration file
3. Modifying existing configuration files easily


----

Run the Configurator
====================

To run the ``Configurator``, the ``monsda.yaml`` conda environment must be installed and activated.

.. code-block::

    conda env create -n MONSDA -f ./MONSDA/envs/monsda.yaml
    conda activate MONSDA


Run the Configurator with

.. code-block::

 monsda_configure

To see further options run the Configurator with the - -help flag

.. code-block::

 monsda_configure --help

----

Main Functions
==============

Create a new project or config-file
-----------------------------------

+---------------------------------------------------------------------------------------------------------------------------+
| ``SELECT WORKFLOWS``                                                                                                      |
| :raw-html:`&vrtri;`                                                                                                       |
| ``MAKE CONDITION-TREE``                                                                                                   |
| :raw-html:`&vrtri;`                                                                                                       |
| ``ASSIGN SAMPLES``                                                                                                        |
| :raw-html:`&vrtri;`                                                                                                       |
| ``SET WORKFLOWS``                                                                                                         |
| :raw-html:`&vrtri;`                                                                                                       |
| ``SET THREADS``                                                                                                           |
+---------------------------------------------------------------------------------------------------------------------------+

To create a new project or a new config-file, the Configurator will take you through all necessary steps.
For creating a project you have to enter a path to establish. Note, that your project folder will grow with your results.
Choose a place with enough memory if necessary.


Modify an existing config-file
------------------------------

1. Add workflows
````````````````

+------------------------------------------------------------------------------------+
| ``SELECT CONFIG``                                                                  |
| :raw-html:`&vrtri;`                                                                |
| ``SELECT ADDITIONAL WORKFLOWS``                                                    |
| :raw-html:`&vrtri;`                                                                |
| ``SET WORKFLOWS``                                                                  |
+------------------------------------------------------------------------------------+

Select Workflows not activated in an existing config-file. The Configurator will
expand each condition. Afterwards you have to make settings for the new workflows.

2. Remove workflows
```````````````````

+------------------------------------------------------------------------------------+
| ``SELECT CONFIG``                                                                  |
| :raw-html:`&vrtri;`                                                                |
| ``SELECT REMOVABLE WORKFLOWS``                                                     |
+------------------------------------------------------------------------------------+

The Configurator will show you all established workflows. After selecting the ones
to be removed it will delete them from the config-file for each condition.

3. Add conditions
`````````````````

+-----------------------------------------------------------------------------------------------------+
| ``SELECT CONFIG``                                                                                   |
| :raw-html:`&vrtri;`                                                                                 |
| ``MAKE CONDITION-TREE``                                                                             |
| :raw-html:`&vrtri;`                                                                                 |
| ``ASSIGN SAMPLES``                                                                                  |
| :raw-html:`&vrtri;`                                                                                 |
| ``SET WORKFLOWS``                                                                                   |
+-----------------------------------------------------------------------------------------------------+

You can add conditions in a similar way you created the condition-tree. But note, that you can't
add sub-conditions to existing leafs. The configurator will expand the condition-tree
for the settings-block and each workflow. Because now you have new option fields in the config-file the Configurator will ask you for copying existing workflow settings or to make new ones.

4. Remove conditions
````````````````````

+-------------------------------------------------------------------------------+
| ``SELECT CONFIG``                                                             |
| :raw-html:`&vrtri;`                                                           |
| ``SELECT REMOVABLE CONDITIONS``                                               |
+-------------------------------------------------------------------------------+

The Configurator will offer you all conditions the condition-tree represents.
After selecting one or several to be removed it will delete them in the
settings-block and for each condition.

----

Interrupt Configuration
=======================

It can happen, that the Configurator asks for entries, you haven't thought about yet.
In this case you can interrupt the configuration and the ``Configurator`` will cache your entries.
A temporary backup file called ``unfinished_config.pkl`` is created for that. 

In most cases you can even just abort the script, but to guarantee clean re-entry you should type

.. code-block::

    exit

When you start the Configurator again later and it finds the ``unfinished_config.pkl`` in the current directory, it will serves a fourth option to continue the session.

Note, that the ``unfinished_config.pkl`` will always be overwritten. To avoid this, you can rename the file.
You can than continue with the --session flag. Run the Configurator like this:

.. code-block:: bash

    monsda_configure -s my_renamed_unfinished_config.pkl

----

Assistance in detail
====================

Create Condition-Tree
---------------------

.. code-block::

  ============================================================

  {
        "NewExperiment": {
              "wild-type": {
                    "day1": {},
                    "day2": {}
              },
              "knockout": {
                    "day1": {},
                    "day2": {}    <=(add sub-conditions here)
              }
        }
  }

  ============================================================

MONSDA understands your experimental design by creating a condition-tree.
The Configurator helps you to create it. To do this, the Configurator points to a condition in which you are allowed to add further sub-conditions.
In this way you can create a nested condition-tree.
Note that each leaf of this tree represents a separate condition. later you will make the workflow settings for each of these conditions.


Sample Assignment:
------------------


.. code-block::

    ============================================================

    {
        "NewExperiment": {
              "wild-type": {
                    "day1": {
                          "SAMPLES": [
                                "Sample_1",
                                "Sample_2"
                          ]
                    },
                    "day2": {}           <-
              },
              "knockout": {
                    "day1": {},
                    "day2": {}
              }
        }
    }

  ============================================================

       1  >  Sample_1     in  path/to/knockout/samples
       2  >  Sample_2     in  path/to/knockout/samples
       3  >  Sample_3     in  path/to/knockout/samples
       4  >  Sample_4     in  path/to/knockout/samples
       5  >  Sample_a     in  path/to/wild-type/samples
       6  >  Sample_b     in  path/to/wild-type/samples
       7  >  Sample_c     in  path/to/wild-type/samples
       8  >  Sample_d     in  path/to/wild-type/samples

The Configurator helps you to assign samples to conditions. If you have activated the FETCH workflow, it will ask you for SRA Accession Numbers.
Otherwise you have to add paths where your samples are stored. The Configurator finds every file with ".fastq.gz" ending and presents it for assignment.
At the same time, the condition-Tree is displayed with an arrow indicating the condition to which samples are assigned.



Make Settings for Conditions
----------------------------

.. code-block::

    ============================================================

      {
            "NewExperiment": {
                  "wild-type": {
                        "day1": {},           <-  1
                        "day2": {},           <-  1
                        "day3": {}            <-  1
                  },
                  "knockout": {
                        "day1": {},           <-    2
                        "day2": {},           <-    2
                        "day3": {}            <-    2
                  }
            }
      }

    ============================================================

MONSDA can run the same workflow with different settings, differentiated by conditions.
Therefore the config-file needs workflow settings for each condition you created.
However you will often set the same settings. To avoid these repetitions during config-creation
the configurator offers you to set several conditions at once.
In the example shown above, you would go through two setting loops.
All sub-conditions of both "wild-type" and "knockout" are assigned the same settings.
To change the conditions set simultaneously, you can loop through the possible selections by pressing enter.
