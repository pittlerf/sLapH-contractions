Author: Martin Ueding

Date: 2018-05-18

---

# Job Script Generator

The contraction code is usually used with a single gauge configuration per
computer job. This has the advantage that there are small work packets that can
be distributed on all the nodes of a cluster. Creating a job script per
configuration is tedious, therefore the contraction code ships with a job
script generator.

The generator is written in Python 3 and uses the Jinja template library for
creating files. It supports `qbig` and `jureca`.

At the time of writing, the help message is the following:

    $ job-script-generator/generate-contraction-jobs -h
    usage: generate-contraction-jobs [-h] [--ignore-checks]
                                     [--conf-skip CONF_SKIP [CONF_SKIP ...]]
                                     --machine {qbig,jureca} --rundir RUNDIR
                                     --outdir OUTDIR --exe EXE [--jobname JOBNAME]
                                     [--jobcount JOBCOUNT] [--email EMAIL]
                                     [--jobs-per-node JOBS_PER_NODE]
                                     [--cpus-per-task CPUS_PER_TASK]
                                     contractinput conf_start conf_end [conf_step]

    Generates a hierarchy of input files and job scripts for the contraction code.
    The script will also make sure that the referenced files actually exist such
    that the jobs will have a higher chance of succeeding.

    positional arguments:
      contractinput         Configuration file contract

    optional arguments:
      -h, --help            show this help message and exit
      --ignore-checks       Do not run the tests for existence of input files.

    Configuration:
      Options that are inherent for the underlying gauge configurations.

      conf_start            First configuration, inclusive
      conf_end              Last configuration, inclusive
      conf_step             default: 1
      --conf-skip CONF_SKIP [CONF_SKIP ...]
                            Skip the given gauge configurations.

    Job:
      Options for the jobs to create.

      --machine {qbig,jureca}
                            Name of the machine the code shall be executed on
      --rundir RUNDIR       Base path for infiles.
      --outdir OUTDIR       Base path for output.
      --exe EXE             Path to the executable. This will be copied into the
                            hierarchy to prevent accidential overwrites.
      --jobname JOBNAME     Name of the submitted job. Default: contraction
      --jobcount JOBCOUNT   Maximal number of jobs submitted to the queue.
                            Default: 60
      --email EMAIL         Email address to send job notifications to. If this is
                            not given, no emails will be send.
      --jobs-per-node JOBS_PER_NODE
                            Number of configurations started in a job. Default: 1
      --cpus-per-task CPUS_PER_TASK
                            Number of cpu cores per task. Default: 4

You have to specify a base configuration file. This file needs to contain all
the options like paths, quarks and operators used. The configuration number
will be passed to the contraction code via the job script and needs not to be
specified. It will be overwritten but for clarity it is best to just omit it in
the configuration file.

Sample configuration files can be found in the directory `infiles` and the one
from the integration test is located at `tests/integration-L4/test_all.ini`.

The job script generator will generate a directory hierarchy with one directory
per gauge configuration. A `start_runs.sh` will also be created to put all jobs
into the queue.
