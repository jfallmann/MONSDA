Integrating new tools/workflows
================================

In case new tools need to be integrated, please refer to similar tools already implemented or contact us in case nothing similar is available yet. Workflows should always be split in subworkflows that follow the same principle as existing subworkflows, ideally making them reusable for other workflows in the future.

Tools should be easy to integrate, all that is needed is to write a tool and if applicable version specific **.smk** or **.nf** file describing input/output and command line calls as well as a fitting **conda environment yaml** file. Once these are available, they should already be usable and configurable via the **config.json** in the specific workflow section.
