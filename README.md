
# Flow


   Flow is an R package that enables local configuration and execution of
   analysis modules on annotated sets of entities (eg pairs, individuals,
   samples). Jobs can be either deployed locally or to a cluster, then
   monitored and managed. Once jobs complete, their outputs can be
   attached back to their respective entities as annotations for easy
   import back into a database or merging with a flat file table.


   Like in the Broad Institute's Firehose platform
   (https://www.broadinstitute.org/cancer/cga/Firehose), a **job** consists of
   a **task** run on an **entity** (e.g. pair, individual, sample). A **task** wraps
   around a module and binds module arguments to names of entity-specific
   annotations or fixed literals which can represent paths (eg a bam file
   path) or values (eg 200). A task also specifies the binding of module
   outputs to output annotations. A **job** is created by applying a task to a
   set of entities, which correspond to keyed table of entity-specific
   annotations (eg bam_file_wgs, seg_file, etc). Once a job completes, one
   or more output annotations (i.e. paths to output files) are attached to
   the respective entity in an output table. See illustration below:
   
   ![Flow Schema](Flow_schema.png)


#   Setting up entities and tasks


   Entities are stored in a keyed R `data.table` of annotations. This table
   can be pulled down from firehose or fiss and imported into R via the
   `data.table` function `fread()`. It can also be obtained via `fiss_get()` in
   db.R or obtained from a data.frame using `as.data.frame()`. The entities
   data.table must have a key (eg pair_id) and that key must have a unique
   value for each entity / row.


   Tasks are configured via an `.task` file. This is a text file whose first
   (non #-commented) line is a path to a firehose module directory (i.e. a
   directory containing a hydrant.deploy file). The subsequent rows are
   tab delimited with 3 or 4 columns, and specify the input and output
   bindings of a task. The first column of every row is 'input' or
   'output'. If the first column is 'input', the second column specifies
   the module argument name that is being bound, the third column
   specifies the annotation, and fourth column has value 'path' or 'value'
   depending on whether the annotation specifies a path or a value. If the
   first column is output, then the second column is the output annotation
   name and the third column is a regexp specifying how to pull the file
   from the module output directory. See the example below.


   An entity table is combined with an `.task` task configuration to create
   an Job, which is a vectorized R object used to run, manage, query, and
   poll the outputs associated with a set of jobs. Instantiation of an Job
   object creates a bunch of subdirectories (by default under `./Flow/`)
   with the task name as sub-directory and entity names as sub-sub
   directories. One can use `cmd()` or `bcmd()` methods to extract shell
   commands for running the jobs locally or on LSF, or the jobs can be
   launched directly from R via the `run()` or `brun()` methods. As jobs are
   executed locally or on LSF, their outputs will be placed into their
   appropriate entity-specific subdirectories (as in firehose), and any
   output annotations that are attached after job completion will refer to
   files residing under their respective task / entity subfolder.


#   Example


   Here's an example of building a dummy module, configuring it to a task,
   and applying that task to a bunch of entities representing tumor normal
   pairs.


   To get started, (install and) load the Flow R package

       install.packages('devtools')
       library(devtools)
       install_github('mskilab/Flow')
       library(Flow)

   We grab the table of entities from a tab delimited file and set the key
   to pair_id. This table comes with the Flow package.

       entities = fread(system.file('extdata', 'entities.txt', package = 'Flow'))
       setkey(entities, pair_id)

   To get things set up we will set up a directory called ~/FlowExample.
   (Make sure you don't have an important directory called ~/FlowExample,
   and if so just replace this path with another path in all the text that
   follows).

       system('mkdir -p ~/FlowExample')
       setwd('~/FlowExample')

   Now we make a module, lets put it in a modules subdirectory of
   FlowExample. (In practice, you will have a static directory containing
   all your modules, which represent reusable code that you will bind to
   multiple task configurations and apply to many different jobs)

       system('mkdir -p ~/FlowExample/modules/dummymodule')

   The module directory contains all the libraries and code associated
   with a given module. It also contains a flow.deploy file containing a
   one-liner with placeholders that will be executed at runtime.


   Let's create the flow.deploy file in this location

       ~/FlowExample/modules/dummymodule/flow.deploy

   In this file, put one line:

       command: <libdir>dummyscript.sh ${ analysis_id } ${ tumor_bam } ${
       normal_bam } ${ error_rate } ${ panel_of_normals } ${ variant_mask }

   All this file needs is a single line prefaced by "command:" (everything
   else is ignored). In this example, that lin contains a call to
   dummyscript.sh and some placeholder variables that are the module
   inputs. The module inputs are specified using the syntax: `${ INPUTNAME
   }` (spaces included).


   The <libdir> in the flow.deploy file is a reserved word that points to
   the directory where the module (and all of its code and files) is
   sitting. At runtime, the code will be run in a job specific output
   directory, and it will need to have a pointer to the module directory
   e.g. to reference scripts and libraries. For example, here it's used to
   specify the path to dummyscript.sh, which we still need to create.


   Let's make the actual dummy script which we store in an executable file


       ~/FlowExample/modules/dummymodule/dummyscript.sh


   In this file we write a simple script that makes a blank vcf and
   outputs a useless report:


       #!/bin/sh
       echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$1.tumor\t$1.normal" > $1.vcf

       echo "analyze bam files $2 and $3 on pair $1 using dummy algo with error rate $4 and panel of normals $5 
       and variant mask $6" >  $1.report.txt


   If you haven't done so make sure the script file is executable


       system('chmod +x ~/FlowExample/modules/dummymodule/dummyscript.sh')


   **(milestone: we've finished making a module!)**


   Next thing to do is to make a task configuration to wrap around this
   module, and store it in a .task file. Again, the .task file binds the
   module arguments to the columns of the entities data.table and/or
   literal arguments, which can be either paths or values. This
   configuration will also be relatively static across many jobs that you
   will end up running on different entities tables.


   First make a tasks directory - this will store task configurations.


       system('mkdir -p ~/FlowExample/tasks')


   To generate a skeleton `.task` configuration for our dummy module, call


       Module('~/FlowExample/modules/dummymodule/')


   This will output to the screen  a skeleton task configuration that you
   can paste into a document and populate with values (i.e. bind
   <INPUT_BINDING> to an actual value and choose path or value) for the
   fourth column.


       #Module dummymodule ("<libdir>dummyscript.sh <tumor_bam> <normal_bam> <error_rate> <panel_of_normals> <variant_mask> <fla...>")
       ~/FlowExample/modules/dummymodule///
       input      tumor_bam         <INPUT_BINDING>       <(path)|(value)>
       input      normal_bam        <INPUT_BINDING>       <(path)|(value)>
       input      error_rate        <INPUT_BINDING>       <(path)|(value)>
       input      panel_of_normals  <INPUT_BINDING>       <(path)|(value)>
       input      variant_mask      <INPUT_BINDING>       <(path)|(value)>
       input      flags             <INPUT_BINDING>       <(path)|(value)>
       output                       <OUTPUT_ANNOTATION>   <OUTPUT_REGEXP>


   Copy and paste this text to a .task file in the tasks directory:
   "`~/FlowExample/tasks/dummy.task`". Now, fill in this task configuration
   by binding the module inputs to table columns or static values (or just
   paste the fully configured task text below).


       #Module dummymodule ("<libdir>dummyscript.sh <tumor_bam> <normal_bam> <error_rate> <panel_of_normals> <variant_mask> <fla...>")
       ~/FlowExample/modules/dummymodule/
       input      tumor_bam         Tumor_clean_bam_file_wgs                        path
       input      normal_bam        Normal_clean_bam_file_wgs                       path
       input      error_rate        "0"                                             value
       input      analysis_id       pair_id                                         value
       input      panel_of_normals  "~/FlowExample/testdata/panel_of_normals.txt"   path
       input      variant_mask      "~/FlowExample/testdata/mask.bed.gz"            path
       output     vcf               .vcf$
       output     quality_metrics   report.txt$
    
   The first non-'#' row in the `.task` file refers to the module directory.
   Note that there are "input" rows and "output" rows, each which are tab
   delimited. The format is pretty flexible, but you need be sure there is
   at least one tab or at least two spaces between each column.


   Note the difference between the input rows (4 columns) vs the output
   rows (3 columns). The third column of each output row specifies a
   regular expression that will be matched against the files outputted in
   each entity's job output directory when it completes. This will bind
   the job-specific output to the annotation specified in the second
   column. Every line beginning with a '#' is ignored.


   Here, we know that the simple script we wrote will output a file ending
   in ".vcf" and "report.txt", and so we configure the task to search the
   job output directory and attach the module outputs matching these
   regular expressions to the annotation fields $vcf and $quality_metrics.


   Using the `Module()` function is one way to generate a .task file.  You
   can also just do it manually, once you understand this simple format.
   If you have trouble getting this task config to run by cutting and
   pasting, it may be a text formatting issue with your terminal app or
   browser corrupting the tab characters.


   To check if the task config parses run:

       Task('~/FlowExample/tasks/dummy.task')

   If this is breaking (sometimes copy / pasting from browser will mess up
   whitespaces) just copy a correctly formatted version of the .task file
   which is distributed with the package and try again.


       system(paste('cp', system.file('extdata', 'dummy.task', package = 'Flow'), '~/FlowExample/tasks/dummy.task'))


   **(milestone: we've finished making a task!)**


   The above procedure (i.e. writing modules, configuring tasks) is something we
   only need to do once or maybe a few times once a module is relatively
   stable.  The stuff that follows we will doing over-and-over again by
   applying many tasks to many datasets.   Back to work:


   Now that we have a  `.task` file, we combine it with the entities table to
   create an `Job` object.   Again, `entities` is a `data.table` and the only
   requirements are that it is (1) it is keyed by a unique entity id and
   (2) it contains all annotation columns that the task is expecting as
   its entity-specific inputs. Now, we're ready to set up jobs for this
   task x entities combination.


       > jobs = Job('~/FlowExample/tasks/dummy.task', entities)
       
           Noting time stamps of inputs
           for tumor_bam (113 paths)
           [1] "tumor_bam"
           for normal_bam (113 paths)
           [1] "normal_bam"
           making output directories under
           /nethome/mimielinski/FlowExample/Flow/dummy
           initializing output annotations
           Dumping out 113 Job.rds files to subdirectories of
           /nethome/mimielinski/FlowExample/Flow
           Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds
           Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds
           
           not ready
           113
           Warning message:
           In .local(.Object, ...) :
           missing annotations resulting causing 113 jobs to be not ready.
           Breakdown of detailed statuses (with # entities with each specific
           status):
           tumor_bam,normal_bam,panel_of_normals,variant_mask not ready(113)


   We've combined the `.task` config with the entities data.table. The
   output 'jobs' is an `Job` object.

       > jobs

           Job on 113 entities (LUAD-CIP-LU-A08-43-T...) with rootdir
           /nethome/mimielinski/FlowExample/Flow from task dummy using module
           dummymodule version
           Job status: not ready (113)


   This vectorized object keeps track of entity specific inputs and
   outputs, stores local, LSF, and SGE / UGER commands for running the
   job, and contains all the task and module information used to build it


   A few things happen as 'jobs' is instantiated: An Flow directory is
   created in the current working directory. This Flow directory contains
   a subdirectory with the task name (Snowman) and 113 subdirectories, one
   for each entity, which will collect task outputs. `Job` is noting the
   time stamps of file path inputs and creating an output data.table that
   will catch job outputs and populate the appropriate entities with
   annotations as specified by the task config.   The location and name of
   this directory can be modified at the time of instantiation (see ?Job).


   Optional arguments to Job include rootdir (changing root directory from
   cwd/Flow), mem (for specifying max memory in GB associated with LSF
   jobs), and queue (for specifying queue on which to run LSF jobs).


   Looking at the output above, we see that the jobs are not ready because
   the input files are missing, i.e. tumor_bam, normal_bam,
   panel_of_normals, and variant_mask. Flow will check inputs and
   determine whether jobs are ready to launch, and has correctly
   determined that these jobs are not ready.


   To get the jobs ready we need to make some "fake data", dummy files
   that sit in the file paths pointed to by columns of the entities table.

       system('mkdir -p ~/FlowExample/testdata')

       sapply(entities$Tumor_clean_bam_file_wgs, function(x) system(paste('touch', x)))
       sapply(entities$Normal_clean_bam_file_wgs, function(x) system(paste('touch', x)))

   We also have to create the variant mask and panel of normal that are
   pointed to as static paths in the task configuration above.

       system('touch ~/FlowExample/testdata/panel_of_normals.txt')
       system('touch ~/FlowExample/testdata/mask.bed.gz')

   Now we can update the jobs object will re-check paths and assess
   readiness to run.

          > update(jobs)
   
          Checking input date stamps
          for tumor_bam (113 files)
          for normal_bam (113 files)
          for error_rate (1 files)
          for analysis_id (113 files)
          for panel_of_normals (1 files)
          for variant_mask (1 files)
          Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds
   
          ready
          113


   Voilà! The jobs are ready so we can run them. Before we do, some basics
   about the Job object:


   It is vectorized, so it has length and can be subsetted. Each element
   corresponds to entity.


          > length(jobs)
          [1] 113
          > jobs[1:10]
          Job on 10 entities with rootdir /cga/meyerson/home/marcin/temp/Flow from task Snowman using module Snowman, version 2015-04-12 10:37:38


   There is a task (Task object) object associated with 'jobs' that you
   can inspect. When you display the task, you see a printout of the .task
   config file.


          > task(jobs)
                 #Module Snowman [Task: Snowman ] ('<libdir>snow.sh <libdir>snowman_150410 run -t <tumor_bam> -n <normal_bam> -e <error_rate> -p ...')
                 /cga/meyerson/home/marcin/svn/CancerGenomeAnalysis/trunk/analysis_pipeline/genepattern/modules/Snowman//
                 input tumor_bam Tumor_clean_bam_file_wgs path
                input normal_bam Normal_clean_bam_file_wgs path
                 input error_rate '0' value
                input cpus '1' value
                 input analysis_id pair_id value
                 input panel_of_normals
                 '/xchip/gistic/Jeremiah/Projects/Lung/lung_snow24_pon.txt.gz' path
                 input indel_mask
                 '/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz' path
                 input flags '--no-r2c-bam' value
                 output snowman_somatic_vcf .*DATECODE.somatic.sv.vcf
                 output snowman_germline_vcf .*DATECODE.germline.sv.vcf
                 output snowman_somatic_indel_vcf .*DATECODE.somatic.indel.vcf
                 output snowman_germline_indel_vcf .*DATECODE.germline.indel.vcf


   Each job has a status (a la firehose) e.g. 'ready', 'complete',
   'failed', etc.


          > status(jobs)[1:3]

                LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
                 'ready'
                 LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
                 'ready'
                 LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE
                 'ready'

   OK let's run some jobs. How about we run the first three entities
   locally (ie on the current machine where R is running).


          > run(jobs[1:3])

                 Starting dummy on entity LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
                 Starting dummy on entity LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
                 Starting dummy on entity LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE


   After launching jobs let's run update() to update the object with the
   latest information. In addition to checking inputs, update will check
   job status, whether a job successfully completed (by analyzing stdout
   and stderr files and polling the directory for relevant outputs).


          > update(jobs)

                 Checking input date stamps
                 for tumor_bam (113 files)
                 for normal_bam (113 files)
                 for error_rate (1 files)
                 for analysis_id (113 files)
                 for panel_of_normals (1 files)
                 for variant_mask (1 files)
                 Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds

                 completed ready
                 3 110

   Great! 3 jobs are finished. We can quickly subset the jobs object using
   a character syntax which searches jobs statuses by regexp. Here we
   subset only the "completed" jobs.


          > jobs['completed']
   
                 Job on 3 entities with rootdir /nethome/mimielinski/FlowExample/Flow from task dummy using module dummymodule version
                 Job status: completed (3)


   This subsetting syntax is especially useful if jobs fail. We can use
   this syntax to debug job failures, restart on a different queue, or
   other change other job specific parameters.


   In this case we have no failed jobs, but we can examine the outputs
   associated with these "completed" jobs. This returns a data.table with
   output annotations as columns. These are file paths attached from the
   respective output directory of each job.


         > outputs(jobs['completed'])

                 pair_id
                 1: LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
                 2: LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
                 3: LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE
                 vcf
                 1: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL.vcf
                 2: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6.vcf
                 3:/nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE.vcf
                 quality_metrics
                 1:/nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL.report.txt
                 2:/nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6.report.txt
                 3:/nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE.report.txt


   We can examine the error stream associated with the first of these
   completed jobs
   
       > more(err(jobs['completed'][1]))

           Command being timed:
           "/nethome/mimielinski/FlowExample/modules/dummymodule///dummyscript.sh LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
            /nethome/mimielinski/FlowExample/testdata//A08-43-4.bam /nethome/mimielinski/FlowExample/testdata//A08-43-1.bam 0
            /nethome/mimielinski/FlowExample/testdata/panel_of_normals.txt /nethome/mimielinski/FlowExample/testdata/mask.bed.gz"
            User time (seconds): 0.00
            System time (seconds): 0.00
            Percent of CPU this job got: 3%
            Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.02
            Average shared text size (kbytes): 0
            Average unshared data size (kbytes): 0
            Average stack size (kbytes): 0
            Average total size (kbytes): 0
            Maximum resident set size (kbytes): 4944
            Average resident set size (kbytes): 0
            Major (requiring I/O) page faults: 0
            Minor (reclaiming a frame) page faults: 345
            Voluntary context switches: 15
            Involuntary context switches: 4
            Swaps: 0
            File system inputs: 8
            File system outputs: 16
            Socket messages sent: 0
            Socket messages received: 0
            Signals delivered: 0
            Page size (bytes): 4096
            Exit status: 0

   We can also query the runtime report to see memory usage and how long
   each completed job took. This outputs another data.table, now with
   runtime information.


       > report(jobs['completed'])

           pair_id
           1: LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
           2: LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
           3: LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE
           jname
           1: dummy.LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
           2: dummy.LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
           3: dummy.LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE
           out.file
           1: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL//dummy.LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL.bsub.out
           2: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6//dummy.LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6.bsub.out
           3: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE//dummy.LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE.bsub.out
           err.file
           1: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL//dummy.LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL.bsub.err
           2: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6//dummy.LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6.bsub.err
           3: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE//dummy.LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE.bsub.err
           exit_flag term_flag started reported
           1: Successfully completed. NA 2016-03-08 18:15:36 2016-03-08 18:15:36
           2: Successfully completed. NA 2016-03-08 18:15:36 2016-03-08 18:15:36
           3: Successfully completed. NA 2016-03-08 18:15:36 2016-03-08 18:15:36
           hours_elapsed max_mem cpu_time success job_type
           1: 0 0.001236 0 TRUE local
           2: 0 0.001240 0 TRUE local
           3: 0 0.001236 0 TRUE local


   Let's run 10 more jobs, this time parallelizing them locally across 5
   cores.

       > run(jobs['ready'][1:10], mc.cores = 5)

           Starting dummy on entity LUAD-CIP-LUAD-AEIUF-TP-NT-SM-1D1NM-SM-1D1K8
           Starting dummy on entity LUAD-CIP-LUAD-D02326-TP-NT-SM-1UVTT-SM-1UVTU
           Starting dummy on entity LUAD-CIP-LUAD-E00934-TP-NT-SM-1UXCW-SM-1UXCX
           Starting dummy on entity LUAD-CIP-LUAD-E01014-TP-NT-SM-1UXD1-SM-1UXD2
           Starting dummy on entity LUAD-CIP-LUAD-E01217-TP-NT-SM-1UXDZ-SM-1UXE1
           Starting dummy on entity LUAD-CIP-LUAD-S00488-TP-NT-SM-18CX6-SM-18CZW
           Starting dummy on entity LUAD-CIP-LUAD-FH5PJ-TP-NT-SM-1D1NY-SM-1D1NZ
           Starting dummy on entity LUAD-CIP-LUAD-E01317-TP-NT-SM-1UXE4-SM-1UXE5
           Starting dummy on entity LUAD-CIP-LUAD-E01278-TP-NT-SM-1UXE2-SM-1UXE3
           Starting dummy on entity LUAD-CIP-LUAD-QY22Z-TP-NT-SM-1DTY7-SM-1DTY8

   Updating the jobs object again, we'll see additional jobs completed.

       > update(jobs)

           Checking input date stamps
           for tumor_bam (113 files)
           for normal_bam (113 files)
           for error_rate (1 files)
           for analysis_id (113 files)
           for panel_of_normals (1 files)
           for variant_mask (1 files)
           Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds

           completed ready
           13 100


   Finally, let's run 5 more of these jobs, but now on the cluster /
   compute farm. On this system, I only have SGE / UGER, so I will use
   `qrun()`. (If I had LSF installed, I would use brun() instead)


       > qrun(jobs[14:18])

           Deploying 7524204 for entity LUAD-CIP-LUAD-S01302-TP-NT-SM-18CXO-SM-18D1F
           Deploying 7524205 for entity LUAD-CIP-LUAD-S01331-TP-NT-SM-18CY2-SM-18D1S
           Deploying 7524206 for entity LUAD-CIP-LUAD-S01341-TP-NT-SM-18CYB-SM-18D22
           Deploying 7524207 for entity LUAD-CIP-LUAD-S01345-TP-NT-SM-18CXP-SM-18D1G
           Deploying 7524208 for entity LUAD-CIP-LUAD-S01346-TP-NT-SM-18CXJ-SM-18D1A


   After a minute I update the jobs again, we'll see additional jobs
   completed.


       > update(jobs)

           Checking input date stamps
           for tumor_bam (113 files)
           for normal_bam (113 files)
           for error_rate (1 files)
           for analysis_id (113 files)
           for panel_of_normals (1 files)
           for variant_mask (1 files)
           Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds

           completed ready
           13 100


   I notice that no additional jobs have completed. I wonder what's wrong
   I can check job status on the farm for these by using `qjobs()`

       > qjobs(jobs[14:18])
   
           pair_id jobid prior name
           1: LUAD-CIP-LUAD-S01302-TP-NT-SM-18CXO-SM-18D1F 7524204 0.00000 dummy.LUAD
           2: LUAD-CIP-LUAD-S01331-TP-NT-SM-18CY2-SM-18D1S 7524205 0.00000 dummy.LUAD
           3: LUAD-CIP-LUAD-S01341-TP-NT-SM-18CYB-SM-18D22 7524206 0.00000 dummy.LUAD
           4: LUAD-CIP-LUAD-S01345-TP-NT-SM-18CXP-SM-18D1G 7524207 0.00000 dummy.LUAD
           5: LUAD-CIP-LUAD-S01346-TP-NT-SM-18CXJ-SM-18D1A 7524208 0.00000 dummy.LUAD
           user state start.sumit.at queue slots taskid
           1: mimielinski qw 03/08/2016 20:49:06 1 NA
           2: mimielinski qw 03/08/2016 20:49:06 1 NA
           3: mimielinski qw 03/08/2016 20:49:06 1 NA
           4: mimielinski qw 03/08/2016 20:49:06 1 NA
           5: mimielinski qw 03/08/2016 20:49:06 1 NA

   Hmm, these are stuck in queue for longer than I expect (not running
   instantly), maybe I need to relaunch with lower memory requirements
   than standard (memory requirement 4)


   First I kill / delete the jobs I just launched.

       > qkill(jobs[14:18])

           mimielinski has deleted job 7524204
           mimielinski has deleted job 7524205
           mimielinski has deleted job 7524206
           mimielinski has deleted job 7524207
           mimielinski has deleted job 7524208

   I set the memory requirement to 1 and rerun

       > mem(jobs) = 1
       > qrun(jobs[14:18])

   Updating the jobs object, I see that 5 more have completed.

       > update(jobs)

           Checking input date stamps
           for tumor_bam (113 files)
           for normal_bam (113 files)
           for error_rate (1 files)
           for analysis_id (113 files)
           for panel_of_normals (1 files)
           for variant_mask (1 files)
           Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds

           completed ready
           18 95


   Now, feeling confident, I push the button to launch all the remaining
   jobs via SGE.   (Note, I don't need to specify indices here, because
   anything that has already completely will not be re-launched)

       > qrun(jobs)

   Updating again, I see everything has completed

       > update(jobs)

          Checking input date stamps
           for tumor_bam (113 files)
           for normal_bam (113 files)
           for error_rate (1 files)
           for analysis_id (113 files)
           for panel_of_normals (1 files)
           for variant_mask (1 files)
           Caching object to /nethome/mimielinski/FlowExample/Flow/dummy.rds

           completed ready
           113 0


   If you want to look under the hood, you'll see that each job is
   associated with a `cmd` (local command), bcmd (LSF command), and `qcmd`
   (UGER / SGE command) which can be accessed using `cmd()`, `bcmd()`, and
   `qcmd()` respectively.   These each return a named vector of shell
   commands. These commands are the result of populating the arguments of
   the module with entity-specific annotation values and static values
   according to the bindings in the task config, as well as adding
   additional shell commands for job tracking. If you want, you can dump
   these commands out into a text file and source it from the shell.


   Jobs are most easily launched directly from the R command line using
   `run(jobs)` or `brun(jobs)`. Alternatively, you can dump the local or bsub
   shell commands to file and run from the shell eg.

       >writeLines(cmd(jobs), 'cmd.sh')
       >writeLines(bcmd(jobs), 'bcmd.sh')
       >writeLines(qcmd(jobs), 'qcmd.sh')

   Once jobs complete, then they will populate their output annotations
   with values after an update() call to the jobs object. You can access a
   data.table of output annotations using the outputs() method eg


       > outputs(jobs)
   
           pair_id
           1: LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
           2: LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
           3: LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE
           4: LUAD-CIP-LUAD-AEIUF-TP-NT-SM-1D1NM-SM-1D1K8
           5: LUAD-CIP-LUAD-D02326-TP-NT-SM-1UVTT-SM-1UVTU
           ---
           109: LUSC-TCGA-77-6843-TP-NB-SM-26XAG-SM-26XAJ
           110: LUSC-TCGA-85-8052-TP-NB-SM-2XLBV-SM-2XLDI
           111: LUSC-TCGA-85-8277-TP-NB-SM-35ASG-SM-35ASJ
           112: LUSC-TCGA-92-8064-TP-NB-SM-2XLDD-SM-2XLCO
           113: LUSC-TCGA-98-8022-TP-NB-SM-2XLCJ-SM-2XLE4

           vcf
           1: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL/LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL.vcf
           2: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6/LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6.vcf
           3: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE/LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE.vcf
           4: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-AEIUF-TP-NT-S
   M-1D1NM-SM-1D1K8/LUAD-CIP-LUAD-AEIUF-TP-NT-
   SM-1D1NM-SM-1D1K8.vcf
           5: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-D02326-TP-NT-
   SM-1UVTT-SM-1UVTU/LUAD-CIP-LUAD-D02326-TP-NT-
   SM-1UVTT-SM-1UVTU.vcf
           ---
           109: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-77-6843-TP-NB-SM-
   26XAG-SM-26XAJ/LUSC-TCGA-77-6843-TP-NB
   SM-26XAG-SM-26XAJ.vcf
           110: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8052-TP-NB-SM-
   2XLBV-SM-2XLDI/LUSC-TCGA-85-8052-TP-NB-
   SM-2XLBV-SM-2XLDI.vcf
           111: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8277-TP-NB-SM-
   35ASG-SM-35ASJ/LUSC-TCGA-85-8277-TP-NB-
   SM-35ASG-SM-35ASJ.vcf
           112: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-92-8064-TP-NB-SM-
   2XLDD-SM-2XLCO/LUSC-TCGA-92-8064-TP-NB-
   SM-2XLDD-SM-2XLCO.vcf
           113: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-98-8022-TP-NB-SM-
   2XLCJ-SM-2XLE4/LUSC-TCGA-98-8022-TP-NB-
   SM-2XLCJ-SM-2XLE4.vcf
           quality_metrics
           1: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM
   -13WXF-SM-13WWL/LUAD-CIP-LU-A08-43-TP-NT-
   SM-13WXF-SM-13WWL.report.txt
           2: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-S
   M-1D1N5-SM-1D1N6/LUAD-CIP-LUAD-2GUGK-TP-NT-\
   SM-1D1N5-SM-1D1N6.report.txt
           3: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-S
   M-1D1ND-SM-1D1NE/LUAD-CIP-LUAD-5V8LT-TP-NT-
   SM-1D1ND-SM-1D1NE.report.txt
           4: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-AEIUF-TP-NT-S
   M-1D1NM-SM-1D1K8/LUAD-CIP-LUAD-AEIUF-TP-NT-\
   SM-1D1NM-SM-1D1K8.report.txt
           5: /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-D02326-TP-NT-
   SM-1UVTT-SM-1UVTU/LUAD-CIP-LUAD-D02326-TP-NT-\
   SM-1UVTT-SM-1UVTU.report.txt
           ---
           109: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-77-6843-TP-NB-SM-
   26XAG-SM-26XAJ/LUSC-TCGA-77-6843-TP-NB
   SM-26XAG-SM-26XAJ.report.txt
           110: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8052-TP-NB-SM-
   2XLBV-SM-2XLDI/LUSC-TCGA-85-8052-TP-NB
   SM-2XLBV-SM-2XLDI.report.txt
           111: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8277-TP-NB-SM-
   35ASG-SM-35ASJ/LUSC-TCGA-85-8277-TP-NB-\
   SM-35ASG-SM-35ASJ.report.txt
           112: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-92-8064-TP-NB-SM-
   2XLDD-SM-2XLCO/LUSC-TCGA-92-8064-TP-NB-
   SM-2XLDD-SM-2XLCO.report.txt
           113: /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-98-8022-TP-NB-SM-
   2XLCJ-SM-2XLE4/LUSC-TCGA-98-8022-TP-NB-
   SM-2XLCJ-SM-2XLE4.report.txt


   Initially, the annotations will be mostly empty, because the jobs have
   not completed. However as jobs complete (as seen above), the values of
   annotation columns will be populated with file paths after running
   `update()` on the jobs object.


   Once jobs complete, the job output can be easily merged with the
   "master table" of entities and updated using a merge() command. This
   merges the tables on their shared key.

       > merge(entities, jobs)

           pair_id
           1: LUAD-CIP-LU-A08-43-TP-NT-SM-13WXF-SM-13WWL
           2: LUAD-CIP-LUAD-2GUGK-TP-NT-SM-1D1N5-SM-1D1N6
           3: LUAD-CIP-LUAD-5V8LT-TP-NT-SM-1D1ND-SM-1D1NE
           4: LUAD-CIP-LUAD-AEIUF-TP-NT-SM-1D1NM-SM-1D1K8
           5: LUAD-CIP-LUAD-D02326-TP-NT-SM-1UVTT-SM-1UVTU
           ---
           109: LUSC-TCGA-77-6843-TP-NB-SM-26XAG-SM-26XAJ
           110: LUSC-TCGA-85-8052-TP-NB-SM-2XLBV-SM-2XLDI
           111: LUSC-TCGA-85-8277-TP-NB-SM-35ASG-SM-35ASJ
           112: LUSC-TCGA-92-8064-TP-NB-SM-2XLDD-SM-2XLCO
           113: LUSC-TCGA-98-8022-TP-NB-SM-2XLCJ-SM-2XLE4
           Tumor_clean_bam_file_wgs
           1: ~/FlowExample/testdata//A08-43-4.bam
           2: ~/FlowExample/testdata//2GUGKAEN.bam
           3: ~/FlowExample/testdata//5V8LTAZB.bam
           4: ~/FlowExample/testdata//AEIUFAVJ.bam
           5: ~/FlowExample/testdata//D02326109.bam
           ---
           109: ~/FlowExample/testdata//TCGA-77-6843-01A-11D-1945-08.bam
           110: ~/FlowExample/testdata//TCGA-85-8052-01A-11D-2244-08.bam
           111: ~/FlowExample/testdata//TCGA-85-8277-01A-11D-2293-08.bam
           112: ~/FlowExample/testdata//TCGA-92-8064-01A-11D-2244-08.bam
           113: ~/FlowExample/testdata//TCGA-98-8022-01A-11D-2244-08.bam
           Normal_clean_bam_file_wgs
           1: ~/FlowExample/testdata//A08-43-1.bam
           2: ~/FlowExample/testdata//2GUGKN18.bam
           3: ~/FlowExample/testdata//5V8LTN9Y.bam
           4: ~/FlowExample/testdata//AEIUFDXJ.bam
           5: ~/FlowExample/testdata//D02326111.bam
           ---
           109: ~/FlowExample/testdata//TCGA-77-6843-10A-01D-1945-08.bam
           110: ~/FlowExample/testdata//TCGA-85-8052-10A-01D-2244-08.bam
           111: ~/FlowExample/testdata//TCGA-85-8277-10A-01D-2293-08.bam
           112: ~/FlowExample/testdata//TCGA-92-8064-10A-01D-2244-08.bam
           113: ~/FlowExample/testdata//TCGA-98-8022-10A-01D-2244-08.bam
           \
           vcf
           1:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM
   -13WXF-SM-13WWL/LUAD-CIP-LU-A08-43-TP-NT-\
   SM-13WXF-SM-13WWL.vcf
           2:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-S
   M-1D1N5-SM-1D1N6/LUAD-CIP-LUAD-2GUGK-TP-NT-\
   SM-1D1N5-SM-1D1N6.vcf
           3:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-S
   M-1D1ND-SM-1D1NE/LUAD-CIP-LUAD-5V8LT-TP-NT-\
   SM-1D1ND-SM-1D1NE.vcf
           4:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-AEIUF-TP-NT-S
   M-1D1NM-SM-1D1K8/LUAD-CIP-LUAD-AEIUF-TP-NT-\
   SM-1D1NM-SM-1D1K8.vcf
           5:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-D02326-TP-NT-
   SM-1UVTT-SM-1UVTU/LUAD-CIP-LUAD-D02326-TP-NT-\
   SM-1UVTT-SM-1UVTU.vcf
           --- \
           109:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-77-6843-TP-NB-SM-
   26XAG-SM-26XAJ/LUSC-TCGA-77-6843-TP-NB-\
   SM-26XAG-SM-26XAJ.vcf
           110:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8052-TP-NB-SM-
   2XLBV-SM-2XLDI/LUSC-TCGA-85-8052-TP-NB-\
   SM-2XLBV-SM-2XLDI.vcf
           111:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8277-TP-NB-SM-
   35ASG-SM-35ASJ/LUSC-TCGA-85-8277-TP-NB-\
   SM-35ASG-SM-35ASJ.vcf
           112:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-92-8064-TP-NB-SM-
   2XLDD-SM-2XLCO/LUSC-TCGA-92-8064-TP-NB-\
   SM-2XLDD-SM-2XLCO.vcf
           113:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-98-8022-TP-NB-SM-
   2XLCJ-SM-2XLE4/LUSC-TCGA-98-8022-TP-NB-\
   SM-2XLCJ-SM-2XLE4.vcf
           \
           quality_metrics
           1:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LU-A08-43-TP-NT-SM
   -13WXF-SM-13WWL/LUAD-CIP-LU-A08-43-TP-NT-\
   SM-13WXF-SM-13WWL.report.txt
           2:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-2GUGK-TP-NT-S
   M-1D1N5-SM-1D1N6/LUAD-CIP-LUAD-2GUGK-TP-NT-\
   SM-1D1N5-SM-1D1N6.report.txt
           3:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-5V8LT-TP-NT-S
   M-1D1ND-SM-1D1NE/LUAD-CIP-LUAD-5V8LT-TP-NT-\
   SM-1D1ND-SM-1D1NE.report.txt
           4:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-AEIUF-TP-NT-S
   M-1D1NM-SM-1D1K8/LUAD-CIP-LUAD-AEIUF-TP-NT-\
   SM-1D1NM-SM-1D1K8.report.txt
           5:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUAD-CIP-LUAD-D02326-TP-NT-
   SM-1UVTT-SM-1UVTU/LUAD-CIP-LUAD-D02326-TP-NT-\
   SM-1UVTT-SM-1UVTU.report.txt
           --- \
          109:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-77-6843-TP-NB-SM-
   26XAG-SM-26XAJ/LUSC-TCGA-77-6843-TP-NB-\
   SM-26XAG-SM-26XAJ.report.txt
           110:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8052-TP-NB-SM-
   2XLBV-SM-2XLDI/LUSC-TCGA-85-8052-TP-NB-\
   SM-2XLBV-SM-2XLDI.report.txt
           111:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-85-8277-TP-NB-SM-
   35ASG-SM-35ASJ/LUSC-TCGA-85-8277-TP-NB-\
   SM-35ASG-SM-35ASJ.report.txt
           112:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-92-8064-TP-NB-SM-
   2XLDD-SM-2XLCO/LUSC-TCGA-92-8064-TP-NB-\
   SM-2XLDD-SM-2XLCO.report.txt
           113:
   /nethome/mimielinski/FlowExample/Flow/dummy/LUSC-TCGA-98-8022-TP-NB-SM-
   2XLCJ-SM-2XLE4/LUSC-TCGA-98-8022-TP-NB-\
   SM-2XLCJ-SM-2XLE4.report.txt


   The resulting entities table has now been updated and may be "ready"
   for downstream pipelines to add additional columns. This table can be
   used to update a central database or saved as a flat file record of
   analytic progress.


# Quick reference / cheat sheet of methods associated with an Job object 'jobs':

| Method                        | Description                                                                                         | 
| ----------------------------- |:---------------------------------------------------------------------------------------------------:| 
| **run(jobs, mc.cores = 1)**       | Run jobs locally from R with optional mc.cores param to run jobs in parallel across mc.cores CPUs.  |
|                               |                                                                                                     |
| **brun(jobs, mc.cores = 1)**      | Submit jobs to LSF from R with optional mc.cores param to submit jobs <br>in parallel across mc.cores CPUs  |
|                               |                                                                                                     |
| **qrun(jobs, mc.cores = 1)**      | Submit jobs to SGE / UGER from R with optional mc.cores param to submit jobs in parallel across mc.cores CPUs  |
|                               |                                                                                                     |
| **update(jobs)**      | Update status of jobs and populate output annotations by scraping / <br>parsing output directories.  |
|                               |  
| **outputs(jobs)**      | Return data.table of output annotations  |
|                               |                        |
| **inputs(jobs)**       | Return data.table of entity-specific input arguments to module.  |
|                               |                        |
| **report(jobs)**       | Return data.table of detailed run info for each job, including max <br>memory usage, run time, job_type (local vs lsf), etc.  |
|                               |  
| **status(jobs)**       | Get short status of jobs (e.g. complete, failed, ready)  |
|                               |      
| **status.info(jobs)**       | Get more detailed status of jobs, including what inputs and outputs <br>might be missing. |
|                               |      
| **jobs[i]**      | Subset jobs using index. 'i' can also be a character, in which case it<br> is interpreted as a regexp against status(jobs). e.g. jobs['complete']<br> will return complete jobs. |
|                               |      
| **cmd(jobs)**       | Get vector of local shell commands associated with jobs.  |
|                               |      
| **bcmd(jobs)**       | Get vector of LSF commands associated with jobs.  |
|                               |      
| **qcmd(jobs)**       | Get vector of SGE / UGER commands associated with jobs.  |
|                               |      
| **cache(jobs)**       | Caches Job object to standard location (<rootdir>/<taskname>.rds)  |
|                               |      
| **key(jobs)**       | Returns the key (name of the column specifying entity id) associated<br> with jobs object.  |
|                               |      
| **cache(jobs)**       | Caches Job object to standard location (<rootdir>/<taskname>.rds)  |
|                               |      
| **ids(jobs)**       | Returns the ids of entities associated with this job object.  |
|                               |    
| **out(jobs)**       | Returns paths to stdout files associated with this job object.  |
|                               |    
| **err(jobs)**       | Returns paths to stderr files associated with this job object.  |
|                               |    
| **dirs(jobs)**       | Returns list of directory contents for each job, takss additional<br> arguments to base function dir()  |
|                               |    
| **purge(jobs)**       | Resets all job output directories to their initial states (deletes all<br> current job outputs, creates fresh directory tree)  |
|                               |    
| **merge(entities, jobs, force = FALSE)**       | Outputs data.table with output annotations from jobs added / updated to<br> 'master' data.table entities. If a given row of entities points to an<br> older path for a given entity-annotation combo, that path will be<br> updated in the output. <br>Otherwise if the entities table contains a newer path than jobs,<br> a warning is thrown and the path is not overwritten <br> in the output (unless force = TRUE).  |
|                               |    


