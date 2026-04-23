# `nf_xpatial`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving nf_xpatial.

We try to manage the required tasks for nf_xpatial using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

## Contribution workflow

If you'd like to write some code for U-BDS/nf_xpatial, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [nf_xpatial issues](https://github.com/U-BDS/nf_xpatial/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf_xpatial repository](https://github.com/U-BDS/nf_xpatial) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline conventions](#pipeline-contribution-conventions)
4. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

### Lint tests

This pipeline enforces the [guidelines set forth by nf-core](https://nf-co.re/developers/guidelines).
If any failures or warnings are encountered, please follow the listed URL for more documentation.

## Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

- On your own fork, make a new branch `patch` based on `upstream/main` or `upstream/master`.
- Fix the bug, and bump version (X.Y.Z+1).
- Open a pull-request from `patch` to `main`/`master` with the changes.

## Pipeline contribution conventions

To make the `nf_xpatial` code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text. 
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test command in `.github/workflow/ci.yml`.

### Default values

Parameters should be initialised / defined with default values within the `params` scope in `nextflow.config`.

Once there, add to `nextflow_schema.json`.

### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/base.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes/steps of the pipeline. 

The process resources can be passed on to the tool dynamically within the process with the `${task.cpus}` and `${task.memory}` variables in the `script:` block.

### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum required version of nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).
