.. role::  raw-html(raw)
    :format: html


==========================================
CONTROLING NextSnakes WITH Configurator.py
==========================================

| NextSnakes operates from a project folder containing a configuration file.
| To enable easy setup and configuration, we host the executable Configurator.
| The Configurator works as an interactive terminal user interface (TUI).
| It supports you in three main things:

1. Creating the complete initial project folder
2. Creating the configuration file
3. Modifying existing configuration files easily

| To learn more about the project folder `click here`_.

.. _projectfolder.rst

| To learn more about structure and function of the configuration file `click here`_.

.. _theconfig.rst

Run the Configurator
====================

To run the Configurator, the ``nextsnakes.yaml`` conda environent must be installed and activated.
To learn more about conda environment management in NextSnakesWrapper `click here`_.

.. code-block::

   conda env create -f ./NextSnakes/envs/nextsnakes
    conda activate NextSnakes


Run the Configurator.py

.. code-block::

 NextSnakes_configure


Main Functions
==============

When the Configurator is started it provides three basic options. To Create a new project, to create a new config-file or to modify an existing config-file.

With the first two options you enter the `Create Mode`_. In both cases you create a configuration file.
The difference with creating a project is that the Configurator asks you for a destination for the project folder and sets soft links for all your input files.

In the following it will be explained in detail, which options the configurator offers you.


Create Mode
------------

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
+===========================================================================================================================+
| To create a new project or a new config-file, the Configurator will take you through all necessary steps.                 |
| For creating a project you have to enter a path to establish. Note, that your project folder will grow with your results. |
| Choose a place with enough memory if necessary.                                                                             |
+---------------------------------------------------------------------------------------------------------------------------+

Modify Mode
-----------



1. Add workflows
````````````````

+----------------------------------------------------------------------------------+
|                                                                ``SELECT CONFIG`` |
|                                                              :raw-html:`&vrtri;` |
|                                                  ``SELECT ADDITIONAL WORKFLOWS`` |
|                                                              :raw-html:`&vrtri;` |
|                                                                ``SET WORKFLOWS`` |
+==================================================================================+
| Select Workflows not activated in an existing config-file. The Configurator will |
| expand it at each condition. Afterwards you have to set the new workflows.       |
+----------------------------------------------------------------------------------+

2. Remove workflows
```````````````````

+-------------------------------------------------------------------------------------+
|                                                                   ``SELECT CONFIG`` |
|                                                                 :raw-html:`&vrtri;` |
|                                                      ``SELECT REMOVABLE WORKFLOWS`` |
+=====================================================================================+
| to be removed it will delete them from the config-file for each condition.          |
| The Configurator will offer you all established workflows. After selecting the ones |
+-------------------------------------------------------------------------------------+

3. Add conditions
`````````````````

+-----------------------------------------------------------------------------------------------------+
|                                                                                   ``SELECT CONFIG`` |
|                                                                                 :raw-html:`&vrtri;` |
|                                                                             ``MAKE CONDITION-TREE`` |
|                                                                                 :raw-html:`&vrtri;` |
|                                                                                  ``ASSIGN SAMPLES`` |
|                                                                                 :raw-html:`&vrtri;` |
|                                                                                   ``SET WORKFLOWS`` |
+=====================================================================================================+
| You can add conditions in a similar way you created the condition-tree. Just add further            |
| subconditions to existing leafes. Afterwards the configurator will expand the condition-tree        |
| for the settings-block and each workflow. Because now you have new option fields in the config-file |
| the Configurator will ask you for copying existing workflow settings or to make new ones.           |
+-----------------------------------------------------------------------------------------------------+

4. Remove conditions
````````````````````

+-------------------------------------------------------------------------------+
|                                                             ``SELECT CONFIG`` |
|                                                           :raw-html:`&vrtri;` |
|                                               ``SELECT REMOVABLE CONDITIONS`` |
+===============================================================================+
| The Configurator will offer you all conditions the condition-tree represents. |
| After selecting the one or several to be removed it will delete them in the   |
| settings-block and for each condition.                                        |
+-------------------------------------------------------------------------------+

Take a Break
============

It can happen, that the Configurator asks for entries, you haven't thought about yet.
So you don't have to abort the creation to start all over again, you can cache your previous entries.
The COnfigurator will safe all your entries in a backup file called ``unfinished_config.pkl``

Whereever you are, type in the terminal:

.. code-block::

    takebreake

Later, to continue the session enter

.. code-block:: bash

    Configurator.py -s unfinished_config.pkl


Operating Assistance in detail
==============================

Create Condition-Tree
---------------------

.. code-block::

  ============================================================

  {
        "NewExperiment": {
              "wildtype": {
                    "day1": {},
                    "day2": {}
              },
              "knockout": {
                    "day1": {},
                    "day2": {}    <=(add subconditions here)
              }
        }
  }

  ============================================================

NextSnakes understands your experimental design by creating a condition-tree.
The Configurator helps you to create it. To do this, the Configurator points to a condition in which you are allowed to add further sub-conditions.
In this way you can create a nested condtion-tree.
Note that each leaf of this tree represents a separate codition. later you can make settings for each of these conditions.


Sample Assignment:
------------------


.. code-block::

    ============================================================

    {
        "DSM1294asdf": {
              "wildtype": {
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
       5  >  Sample_a     in  path/to/wildtype/samples
       6  >  Sample_b     in  path/to/wildtype/samples
       7  >  Sample_c     in  path/to/wildtype/samples
       8  >  Sample_d     in  path/to/wildtype/samples

enter all sample-numbers according to the displayed condition comma separated
>>> 3,4



Make Settings for Conditions
----------------------------

.. code-block::

    ============================================================

      {
            "NewExperiment": {
                  "wildtype": {
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

In the following steps you will make different settings for each condition.
To avoid repetitions, specify which conditions should get the same settings
You will set all conditions with the same number at once afterwards

To loop through the possible selections press enter
Finally enter 'ok' to make the settings
