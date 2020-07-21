=======================================
CREATE config.json WITH Configurator.py
=======================================

To enable users easy ``snakemake`` configuration, we host the
executable ``Configurator.py`` which generates a ``config.json``
file based on the users needs. This ``config.json`` follows standard
``json`` format and contains all the information needed to run the
analysis.  Making use of ``Configurator.py`` allows to start with a
simple ``config.json`` based on the template in the configs
directory and append workflows as needed step by step, or create the
full configuration at once.  Usually a user will start an analysis
running initial Quality Control on raw files. Although
``RunSnakemake.py`` will run QC (when enabled) also for trimming and
mapping steps, this will be used as first step explaining
``Configurator.py``.

To list all options and available choices for each option simply run:
``
./snakes/Configurator.py
``

To create an initial ``config.json`` simply run:
``
./snakes/Configurator.py -w QC -i hs:01012020:std -c config_basic.json -m hs:hg38 -g hg38:hg38genomefile -x hg38:all_chromosomes
``

Where ``-w`` is the option to add 'QC' as workflow step, ``-i``
sets the IdentifierConditionSetting relationship ``-c`` names the
output file where the configuration is written to and so on, please
refer to the help output of the tool for more information.

After 'QC', the user wants to start trimming and mapping and also add
new ICSs so we append to the just created config to make
``RunSnakemake.py`` aware of that.

``
./snakes/Configurator.py -w QC,TRIMMING,MAPPING -i hs:01012020:std,newlab:12012020:specific,otherspecies:24032003:test -c config_basic.json -m hs:hg38,newlab:hg38,otherspecies:dm6 -g hg38:hg38genomefile,dm6:dm6genomefile -x hg38:all_chromosomes
``

And ``Configurator.py`` will append to the existing
``config.json`` everything that is new and overwrite what has
changed and keep the rest intact.  Following this procedure the user
can stepwise analyze samples one by one, step by step or run the whole
pipeline at once with all samples and processing steps. It is also
possible to append to another existing ``config.json``, simply
change the filename of the ``-c`` option.

When run with the standard template.json as input, most setting will
be very basic, only a subset of tools will be used so we recommend to
manually go over the configuration and adapt settings to your needs.
