##############################################################################
## Marcin Imielinski
##
## Weill Cornell Medicine
## mai9037@med.cornell.edu
##
## New York Genome Center
## mimielinski@nygenome.org
##

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

#' @import parallel
#' @import data.table
#' @importMethodsFrom data.table key
#' @import stringr


#' @name Module-class
#' @title Module class
#' @description
#'
#' Class \code{Module} defines an unconfigured module that is scraped from a "module" or "hydrant.deploy" file.
#'
#' @section Slots:
#' \describe{
#'   \item{name}{character name of module (by default the name of the module directory)}
#'   \item{sourcedir}{character path of firehose module directory containing code and hydrant.deploy file}
#'   \item{args}{vector of argument names that need to be populated for this module to be configured into a task}
#'   \item{cmd}{formatted string to be populated by args at time of instantiation of this module into Job}
#'   \item{stamp}{version stamp of this module}
#' }
#'
#' @exportClass Module
#' @author Marcin Imielinski
setClass('Module', representation(sourcedir = 'character', name = "character", cmd = 'character', args = 'vector', stamp = 'character'))

setMethod('initialize', 'Module', function(.Object,
                                             path, # path to module
                                             name = NULL ## default is 
                                             )
    {
        if (!file.exists(path)) 
            stop(sprintf('%s not found, check path', path))

        if (file.info(path)$isdir)
            path = paste(path, 'hydrant.deploy', sep = '/')
        else
            path = paste(file.dir(path), 'hydrant.deploy', sep = '/')

        if (!file.exists(path))
            path = paste(file.dir(path), 'flow.deploy', sep = '/')
        
        if (!file.exists(path)) 
            stop(sprintf('%s not found, check path', path))

        cmd.re = '^command\\s*[\\:\\=]\\s+'
        cmd = gsub(cmd.re, '', grep(cmd.re, readLines(path), value = TRUE))

        if (length(cmd)==0)
            stop('Problem parsing .deploy file, maybe it is an old version or a scatter gather task, which is not currently supported')

        if (!grepl('sh', cmd))
            warning('Module command does not begin with "sh" - make sure that this is not a scatter gather command, which is not currently supported')

        
        ## need to replace $(\\w+ .* FEATURE_NAME) with just the internal and extract the FEATURE_NAME
        pattern = '\\$\\{[a-z\\,]*( [^\\}]*)? (\\S+)\\s*\\}';              
        args = stringi::stri_match_all_regex(cmd, pattern, omit_no_match = TRUE, cg_missing = "")[[1]]

        args[,3] = gsub('\\=.*', '', args[,3])
        
        ## create new string with clear "formats" we can sub into
        cmd.new = cmd
        for (i in 1:nrow(args))
            cmd.new = str_replace_all(cmd.new, stringr::fixed(args[,1][i]), paste(gsub('\\"', '', args[,2]), '<', args[,3], '>', sep = '')[i])
        
        .Object@cmd = cmd.new

        .Object@args = args[,3]

        .Object@stamp = ''
        if (length(f <- dir(file.dir(path), '^version.txt$'))>0)
            .Object@stamp = readLines(f[1])
        ## else ## extract from svn 
        ##     tryCatch(
        ##         {
        ##             p = pipe(sprintf('svn info %s | grep "Last Changed Date"', normalizePath(file.dir(path))))
        ##             .Object@stamp = str_match(readLines(p)[1], 'Last Changed Date\\: (\\S+\\s+\\S+)')[,2]
        ##             close(p)
        ##         }, error = function(e) warning('Could not get stamp'))        
        .Object@sourcedir = file.dir(path)
        if (is.null(name))
            .Object@name = file.name(gsub('\\/+$', '', file.dir(path)))
        else
            .Object@name = name
            
        return(.Object)
    })

#' @name Module
#' @title Class to implement modules, which are standalone pieces of code that expect standard inputs and outputs
#' @description
#'
#' Module constructor takes as input the path to a directory containing either a .module or hydrant.deploy file (for compatibility with firehose)
#' The module file is scraped to look for a single line beginning with "command:" that specifies the code to be run which placeholders for inputs
#' and specifying code that only writes data to <relative paths> at or below the current working directory.
#'
#' Module('path.to.module.directory.with.module.or.hydrant.deploy.file')
#' 
#'
#' The input placeholders are specified using the following syntax: ${t Argument_Name}.  These inputs will eventually be attached to entity
#' annotation. 
#'
#' There is also a line <libdir> that specifies the file path of the <module directory>.  Since the code will be eventually executed in a
#' Job specific output (hence the need for the module to only write to relative paths) there needs to be a way
#' the code to refer to other files in the module directory.  This is what <libdir> provides a handle to. 
#'
#' So an example module file for JabbA has a single line:
#' command: sh <libdir>run.sh <libdir>run.jabba.R -l <libdir> -n ${t TumorName} -s ${t SegFile} -a ${t CovFile} -r ${t RAfile} -g ${t SubSample} -b ${t NormalSegFile} -k ${t SlackPenaltyPerLooseEndCopy} --iterate ${t NumIterations} --tfield ${t TierFieldName} --hets ${t OptionalHetPileupOutput} 
#'
#' 
#' Now, in most cases run.sh will be a pretty generic top level shell script that will set up the environment (eg loading additional libraries, setting environment variables) and pipe the remaining commands to an R, python, or perl script.  In other cases this can be a script doing more of the "heavy lifting" in the task.  This is up to the user.  However the .module file itself should be a one-liner as above. 
#'
#' Regarding outputs, The module definitoin does not currently specify what are its output files (this may change).
#' These are currently specified at the task level, where the task author (who presumably knows the module) will specify regexp to scrape
#' specific files out of the Job output directory and attach their paths to the respective entity. 
#' 
#' @param path character path to module directory containing .module or .deploy file
#' @author Marcin Imielinski
#' @export
Module = function(...) new('Module', ...)

setMethod('show', 'Module', function(object)
    {
        writeLines(as.character(object))
    })

setMethod('as.character', 'Module', function(x, ...)
    {
        NCHAR = 100
        tmp.cmd = paste(str_sub(x@cmd, 1, pmin(nchar(x@cmd), NCHAR)), ifelse(nchar(x@cmd) > NCHAR, '...', ''), sep = '')
#        cat(sprintf('module %s, with command:  %s\nand args:\n%s\n', x@name, tmp.cmd, paste(x@args, collapse = '\n')))
        out = sprintf('#Module %s ("%s")', x@name, tmp.cmd)
        out = c(out, x@sourcedir, paste('input', x@args, '<INPUT_BINDING>', '<(path)|(value)>', sep = '\t'))
        out = c(out, paste('output\t<OUTPUT_ANNOTATION>\t<OUTPUT_REGEXP>'))
        return(out)
    })


#' S4 class for \code{Task}
#'
#' Class \code{Task} defines an task configuration object that is built from a firehose module with a 
#' hydrant.deploy file, and can be run locally or on the farm on annotations provided in a
#' data.table 
#'
#' @section Slots:
#' \describe{
#'   \item{name}{scalar character task name}
#'   \item{module}{Module object}
#'   \item{libdir}{scalar character directory to the FH libdir (by default is the module directory)}
#'   \item{outdir}{scalar character corresponding to name of out directory}
#'   \item{mem}{numeric scalar specifying default maximum memory in G to allocate to this task when runnin on farm}
#'   \item{args}{list of arguments, each is either an FlowLiteral or FlowAnnotation this specifies how the Task will be instantiated into an Job when combined with a data.table, can also be a character in which is interpreted as an annotation name by default}
#'   \item{outputs}{list of FlowOutput objects}
#'   \item{stamp}{timestamp of this task configuration instantiation}
#'  }
#'
#' @name Task-class
#' @rdname Task-class
#' @exportClass Task
#' @export
setClass('Task', representation(name = 'character', module = 'Module', mem = 'numeric', libdir = 'character', args = 'list', stamp = "character", outputs = "list"))


#' @export
setMethod('initialize', 'Task', function(.Object,
                                           config, ## can be Module or character path to text based Task config file
                                           mem = 4,
                                           name = NULL,
                                           libdir = NULL,
                                           output = NULL, ## FlowOutput object or list of FlowOutput objects
                                           ... ## arguments to the module, and their values are either FlowLiteral or FlowAnnotation.  Can also be a scalar character in which case interpreted as an annotation name
                                           )
    {
        module = config

        if (is.character(module)) ## interpret as path
            {
                errstr = 'First non-#-commented line of task config file must be path to module directory with .deploy file, and remaining lines must consist of three or more tab delimited columns.  \nThe first column is in each line specifies whether that line is either input or output (with the corresponding text),  \nThe second column is a module argument name (if input) or an output annotation name (if output).  \nIf the line is an input, the third column is an unquoted string specifying an annotation name or a quoted string specifying a literal.  \nIf the line is an output, the third column is a bare string specifying a regex which tells Flow where to find the file to attach to the annotation\nFor inputs, an optional fourth column can specify whether it is a "path" or "value".  For input annotations an optional 5th column can specify default values'
                if (!file.exists(module))
                    stop(paste('module file', module, 'does not exist'))

                lines = strsplit(grep('\\S+',
                    grep('^#', readLines(module), invert = TRUE, value = TRUE),
                    value = TRUE), '(\\s*\t\\s*)|(\\s\\s+)')
                lens = unlist(sapply(lines, length))
                if (lens[1] != 1)
                    stop(paste('Error reading .task task config file:\n', errstr))

                if (length(lens)>1)
                    if (!all(lens[-1]>=3))
                        stop(paste('Error reading module file:\n', errstr))
                
                modfn = str_trim(lines[[1]])
                lines = lines[-1]

                if (length(lines)>0)
                    {
                        tmp = sapply(lines, function(x) x[1])
                        ltype =ifelse(grepl('input', tmp, ignore.case = TRUE), 'input',
                            ifelse(grepl('output', tmp, ignore.case = TRUE), 'output', 'other'))
                        
                        if (any(ltype=='other'))
                            stop(paste('Some input / output lines have invalid first column:\n', errstr))
                    }
                
                if (!file.exists(modfn))
                    stop(paste('Module path', modfn,  'pointed to by this task config file does not exist.  Check the path and read format spec belwow:\n', errstr))
                
                .Object@module  = Module(modfn)
                .Object@libdir  = .Object@module@sourcedir
                .Object@mem  = as.numeric(mem)

                name = gsub('\\.task$', '', file.name(module))
                
                tmp.args = NULL
                output = NULL
                
                if (length(lines)>0)
                    {
                        lqre = '^\\s*[\\"\\\']'
                        rqre = '[\\"\\\']\\s*$'
                        .gsub = function(x) gsub(rqre, '', gsub(lqre, '', x))
                        if (any(ltype == 'input'))
                            {
                                tmp.args = lapply(lines[ltype=='input'], function(x, nm)
                                    {
                                        if (grepl(lqre, x[3]) | grepl(rqre, x[3]))
                                            {
                                                path = ifelse(is.na(x[4]), FALSE,
                                                    ifelse(grepl('(true)|(path)', x[4], ignore.case = TRUE), TRUE, FALSE))
                                                FlowLiteral(x[2], .gsub(x[3]), path)
                                            }
                                        else if (suppressWarnings(!is.na(as.numeric(x[[3]]))))
                                            {
                                                FlowLiteral(x[2], x[3], FALSE)
                                            }
                                        else
                                            {
                                                path = ifelse(is.na(x[4]), TRUE,
                                                    ifelse(grepl('(false)|(val)', x[4], ignore.case = TRUE), FALSE, TRUE))

                                                default = NA
                                                if (!is.na(x[5]))
                                                    default = gsub(rqre, '', gsub(lqre, '', x[5]))
                                                FlowAnnotation(x[2], x[[3]], path, default)
                                            }
                                    })
                            }
                        
                        names(tmp.args) = sapply(lines[ltype=='input'], function(x) x[2])
                        
                        if (any(ltype == 'output'))
                           output = lapply(lines[ltype=='output'], function(x) FlowOutput(x[2], .gsub(x[3])))
                    }
            }
        else            {                
                .Object@module = module
                
                if (is.null(libdir))
                    libdir = module@sourcedir
                
                .Object@libdir = libdir        
                .Object@mem = as.numeric(mem)
        
                if (is.na(.Object@mem))
                    stop('Must provide numeric value for default mem')
                
                tmp.args = list(...)
            }

        if (is.null(name))
            name = .Object@module@name
        
        .Object@args = lapply(tmp.args, function(arg)
            {
                if (arg@name == 'job.spec.memory' & is(arg, 'FlowLiteral'))
                    {
                        is.num = suppressWarnings(!is.na(as.numeric(arg@arg)))
                        if (is.num)
                            .Object@mem <<- as.numeric(arg@arg)                        
                    }
                else if (!(arg@name %in% .Object@module@args))
                    warning(paste(arg@name, 'is not an argument to', .Object@module@name))

                
                if (is(arg, 'character'))
                    {
                        if (length(arg)>1)
                            stop('argument must be length 1 character vector or FlowAnnotation / FlowLiteral objects')
                        arg = FlowAnnotation(arg)
                    }


                return(arg)
            })
        names(.Object@args) = names(tmp.args)

        if (length(missing <- setdiff(.Object@module@args, names(.Object@args)))>0)
            stop(sprintf('These module arguments are not specified in the provided task configuration:\n%s\nplease fix ...', paste(missing, collapse = ', ')))

          if (!is.list(output))
              output = list(output)

          if (!all(sapply(output, is.null) | sapply(output, is, 'FlowOutput')))
              stop('output arg must be FlowOutput, list of FlowOutput objects, or NULL')
          
          .Object@outputs = output
          .Object@name = name


          return(.Object)
     })



#' @name Task
#' @title Object representing a task, which is a wiring of a module inputs and outputs to specific entity annotations names
#' @description
#'
#' initializing Task
#' requires Module object (or path to module file), and values for all arguments required by Task with literals or names of columns
#' from which to populate eventual job data. 
#' also can specify default mem for task (can be changed at job level)
#' output is a task configuration
#' module can also be a text based task configuration file which will have the format:
#' MODULE_PATH (
#' input       [(annotation)|(literal)]         VALUE
#' ...
#' ...
#' output      ANNOTATION_NAME                  REGEXP_FOR_FILE
#'
#' ie first row is the module path, the next rows specify whether input or output is being specified
#' in the first column, and if input, the second column specifies whether annotation or literal and third column specifis value
#' otherwise for output, second column specifies the name of the output annotation and the third column specifies the
#' regexp used to query the output directory for the output file 
#' 
#' a task object is normally instantiated from a text .task file.
#'
#' It wraps together Module, FlowLiteral, FlowOutput, and FlowAnnotation objects and can also be instantiated directly from them.
#'
#' The task definition file has a header, pointer to a directory containing all the code necessary to run the module + a .deploy
#' file (a la firehose .deploy), and inputs / outputs in a standard syntax.
#'
#' Here is an example of a task  definition for Jeremiah Wala's Snowman rearrangement detector
#' 
#' #Module Snowman ('<libdir>snow.sh <libdir>snowman_150410 run  -t  <tumor_bam>  -n  <normal_bam>  -e  <error_rate>  -p ...')
#' ~/modules/Snowman/
#' input   tumor_bam       Tumor_clean_bam_file_wgs        path
#' input   normal_bam      Normal_clean_bam_file_wgs       path
#' input   error_rate      '0'     value
#' input   cpus    '1'     value
#' input   analysis_id     pair_id value
#' input   panel_of_normals        '/xchip/gistic/Jeremiah/Projects/Lung/lung_snow24_pon.txt.gz'   path
#' input   indel_mask      '/xchip/gistic/Jeremiah/Projects/HengLiMask/um75-hs37d5.bed.gz'         path
#' input   flags   '--no-r2c-bam'                                                                  value
#' output  snowman_somatic_vcf     .*DATECODE.somatic.sv.vcf
#' output  snowman_germline_vcf    .*DATECODE.germline.sv.vcf
#' output  snowman_somatic_indel_vcf       .*DATECODE.somatic.indel.vcf
#' output  snowman_germline_indel_vcf      .*DATECODE.germline.indel.vcf
#'
#' 
#' @param config path to text based ".task" task config file or Module object 
#' @param mem memory limit to task (default 4)
#' @param name name of task
#' @param libdir libdir library / module directory
#' @param output  FlowOutput object or list of FlowOutput object
#' @param ...  additional FlowLiteral or Flow Annotation boject
#' @author Marcin Imielinski
#' @export
Task = function(...) new('Task', ...)

setMethod('show', 'Task', function(object)
    {
        writeLines(as.character(object))
    })

#' @name as.character
#' @title Convert task object to character, can be dumped to file and used to instantiate future Task objects
#' @export
#' @author Marcin Imielinski
setMethod('as.character', 'Task', function(x, ...)
    {
        NCHAR = 100
        tmp.cmd = paste(str_sub(x@module@cmd, 1, pmin(nchar(x@module@cmd), NCHAR)), ifelse(nchar(x@module@cmd) > NCHAR, '...', ''), sep = '')        
        out = sprintf('#Module %s [Task: %s ] ("%s")', x@module@name, x@name, tmp.cmd)
        if (length(x@args)>0)
            out = c(out, x@module@sourcedir,
                paste('input\t',
                      names(x@args),'\t', 
                      ifelse(sapply(x@args, is, 'FlowLiteral'), '"', ''),
                      sapply(x@args, function(x) x@arg),
                      ifelse(sapply(x@args, is, 'FlowLiteral'), '"', ''),
                      '\t',
                      ifelse(sapply(x@args, function(x) x@path), 'path', 'value'),
                      '\t',
                      ifelse(sapply(sapply(x@args, default), is.null), '', sapply(x@args, default)),
                      sep = ''))
        
        if (length(x@outputs)>0)
            out = c(out,
                paste('output\t',
                      sapply(x@outputs, function(x) x@name), '\t',
                      sapply(x@outputs, function(x) x@pattern)))

        return(out)
    })


#' @name FlowLiteral-class
#' @title Class to represent a literal argument to a task (ie a static path or value), used for wiring modules into tasks 
#' @description
#' S4 class for \code{FlowLiteral}
#'
#' Class \code{FlowLiteral} is a simple container for storing arguments to Task that are to be interpreted literally
#' used when initializing a task config.
#'
#' @section Slots:
#' \describe{
#'   \item{arg}{scalar character argument that will be interpreted literally}
#' }
#' 
#' @rdname FlowLiteral-class
#' @exportClass FlowLiteral
#' @export
#' @author Marcin Imielinski
setClass('FlowLiteral', representation(name = 'character', arg = 'character', path = 'logical'))
setMethod('initialize', 'FlowLiteral', function(.Object,
                                              name, # argument name
                                              arg,  # argument value
                                              path = FALSE # is it a path or not
                                              )
    {
        .Object@name =  name
        .Object@arg =  arg
        .Object@path =  path
        return(.Object)
    })



#' @name FlowLiteral
#' @title Creates object representing a literal argument to a task (ie a static path or value), used for wiring modules into tasks 
#' @param name character name of argument
#' @param arg value of argument
#' @param path boolean specifyign whether this is a path or not
#' @export
#' @author Marcin Imielinski
FlowLiteral = function(...) new('FlowLiteral', ...) 

#' @name FlowAnnotation-class
#' @title Class to represent an entity-specific annotation (ie a path or value), used for wiring modules into tasks 
#' @description
#' 
#' Class \code{FlowAnnotation} is a simple container for storing arguments to Task that are to be interpreted as annotations
#' used when initializing a task config.
#' S4 class for \code{FlowAnnotation}
#'

#' @section Slots:
#' \describe{
#'   \item{arg}{scalar character argument that will be interpreted literally}
#' }
#' @exportClass FlowAnnotation
#' @export
setClass('FlowAnnotation', representation(name = 'character', arg = 'character', path = 'logical', default = 'character'))

setMethod('initialize', 'FlowAnnotation', function(.Object,
                                                 name, # argument name 
                                                 arg, # argument value,
                                                 path = TRUE, # flag whether it is path
                                                 default = NULL ## specify default value
                                             )
    {
        .Object@name =  name
        .Object@arg =  arg
        .Object@path =  path
        
        if (is.na(default))
            default = as.character(NULL)
        
        .Object@default = default
        return(.Object)
    })


#' @name FlowAnnotation
#' @title Creates an entity-specific annotation object (ie a path or value), used for wiring modules into tasks 
#' @param name character name of argument
#' @param arg value of argument
#' @param path boolean specifying whether this is a path or not
#' @param default can specify default value (ie what do to if annotation is missing from input table) (default is NULL)
#' @export
FlowAnnotation = function(...) new('FlowAnnotation', ...) 


#' @name FlowOutput-class
#' @title Class representing the wiring of a regexp on a job output directory to an output annotation

#' S4 class for \code{FlowOutput}
#'
#' Class \code{FlowOutput}is a simple container for capturing  output from Tasks
#' which consists of name (of output annotation) and pattern (to match the output path in the
#' task's output directory)
#'
#' @section Slots:
#' \describe{
#'   \item{name}{scalar character specifying annotation to write}
#'   \item{pattern}{scalar character regexp pattern to match the file in the output directory}
#' }
#' @name FlowOutput-class
#' @rdname FlowOutput-class
#' @exportClass FlowOutput
#' @export
setClass('FlowOutput', representation(name = 'character', pattern = 'character'))
setMethod('initialize', 'FlowOutput', function(.Object,
                                             name,
                                             pattern
                                             )
    {
        .Object@name =  name
        .Object@pattern = pattern
        return(.Object)
    })

#' @name FlowOutput
#' @title  Constructs an object representing the wiring of a regexp on a job output directory to an output annotation
#' @description
#' 
#' Constructs an object representing the wiring of a regexp on a job output directory to an output annotation
#' 
#' @param name annotation name to attach output file
#' @param pattern pattern to match to get output file in output directory
#' @export
#' @author Marcin Imielinski
FlowOutput = function(...) new('FlowOutput', ...) 


#' @name Job-class
#' @rdname Job-class
#' @title Class representing a set of jobs, each which is the combination of a Task and entity (i.e. a sample, pair, or set)
#' @description
#' S4 class for \code{Job}
#'
#' Class \code{Job} defines a set of nrow(entities) firehose jobs "myjobs" which is an instantiation of an Task combined with
#' one or more rows of an (keyed) annotation data.table entities.  An Job is created from an Task + data.table of annotations, and (optional)
#' output directory (otherwise made default).   Once created, several things can be done with an Job:  It's local command line tasks
#' can be extracted with cmd(myjobs), it's bsub tasks can be extracted with bcmd(myjobs). These can be dumped to a file or executed 
#' directory from R via a system command on the cmd(), bcmd(), qcmd() output, or run with run(myjobs),  brun(myjobs), qrun(myjobs)...  
#'
#' Once jobs are running / completing the status of jobs can be checked with update(.Object).
#' In this case the timestamp of each output will be updated and the relevant annotations will be
#' updated after querying the output directories.  The data.table associated with the Job can be accessed
#' using the runinfo command i.e. runinfo(.Object)
#'
#'
#' @section Slots:
#' \describe{
#'   \item{task}{Task object that this job was built from}
#'   \item{rootdir}{scalar character directory specifying root output directory of jobs}
#'   \item{runinfo}{keyed data.table with length(.Object) rows containing columns
#'            $outdir = output directory for each job
#'            $queue = character specifying the queue on which this job should be run
#'            $mem = numeric specifying G of max memory to run this job with 
#'            $cmd = command line to execute job locally
#'            $bcmd = command line to submit job to farm
#'            $status = character that is either "completed", "ready", "not ready", or "outdated" if time
#'            stamp of one or more input annotation >= timestamp of the output annotation.
#'            $status.info = description of status when not ready, which can describe missing inputs or updated inputs since last run}
#'   \item{inputs}{keyed data.table with length(.Object) rows with column corresponding to key and columns specifying values
#'             input annotations for each entity}
#'   \item{stamps}{keyed data.table with length(.Object) rows with column corresponding to key and columns specifying timestamps
#'             of input annotations for each entity}
#'   \item{outputs}{keyed data.table with length(.Object) rows $timestamp column for time stamp of output (or NA if job is not
#'            complete) and columns for table key + 0 or more output annotations with filled in values if job is complete or NA}
#' }
#' @exportClass Job
#' @export
setClass('Job', representation(task = 'Task', rootdir = 'character', entities = 'data.table', runinfo = 'data.table', inputs = "data.table", stamps = "data.table", outputs = "data.table"))

setMethod('initialize', 'Job', function(.Object,
                                          task, ##Task wrapping around an Module expecting literal and annotation arguments
                                                ##this can also just be an .task config file
                                          entities = NULL, ## keyed data.table, key will determine id of outgoing jobs, columns of table used to populate task
                                          rootdir = './Flow/',
                                          queue = as.character(NA),
                                          nice = NULL,
                                          mem = NULL,                                          
                                        cores = 1,
                                        now = FALSE,
                                        mock = FALSE
                                          )
    {
        require(stringr)
        if (is.null(nice))
            nice = TRUE

        if (is.character(task)) ##
            task = Task(task)
        else if (!is(task, 'Task'))
            stop('Job must be instantiated from Task object or .fhtask file')

        if (!is.data.table(entities))
            {
                if (is.data.frame(entities))
                    {
                        id = rownames(entities)
                        if (is.null(id))
                            stop('Input must be keyed data.table or data.frame with rownames')                         
                        warning('Converting entities data.frame to data.table setting key "id" from rownames')
                        entities = as.data.table(entities)
                        entities$id = id
                        setkey(entities, 'id')                                          
                    }
                else
                    stop("Entities input must be keyed data.table or data.frame with non NULL and unique rownames'")
            }
        else
            if (is.null(data.table::key(entities)))
                stop('Entities argument must be a keyed data.table, please add a key')

        if (any(duplicated(entities[[data.table::key(entities)]])))
            stop(sprintf('Input entities table has duplicated keys! Check column "%s" of entities table or use a different column as a key', data.table::key(entities)))
        
        .Object@entities = copy(entities)

        tabstring = function(tab, sep = ', ')
            return(paste(names(tab), '(', tab, ')', sep = '', collapse = sep))        

        ann.args = do.call(c, lapply(task@args, function(x)
            if (is(x, 'FlowAnnotation'))
                return(x)
            else
                return(NULL)))

         tmp = task@args[sapply(task@args, is, 'FlowLiteral')]
        lit.args = structure(lapply(tmp, function(x) x@arg), names = names(tmp))

        if (any(iz <- sapply(ann.args, function(x) !is.null(default(x)))))
            if (any(ix <- !sapply(ann.args[iz], function(x) x@arg) %in% names(entities)))
                for (arg in ann.args[iz][ix])
                    entities[[arg@arg]] = NA        
        
         if (!all(ix <- sapply(ann.args, function(x) x@arg) %in% names(entities)))
             stop(sprintf('These annotations / columns required for task %s are missing from entities data.table: %s', task@name, paste(sapply(ann.args, function(x) x@arg)[!ix], collapse = ', ')))

         if (any(duplicated(entities[, data.table::key(entities), with = FALSE][,1])))
             stop('Duplicate entities present with respect to current key.  Each entity should have a unique key.  Check entities table and/or key settings')
         
         module = task@module        
         .Object@task = task

        if (!mock)
            {
                system(paste('mkdir -p', rootdir))
                rootdir = normalizePath(rootdir)                                
            }
        .Object@rootdir = rootdir
        .Object@stamps = data.table(entities[, data.table::key(entities), with = FALSE])
        .Object@outputs = data.table(entities[, data.table::key(entities), with = FALSE])
        .Object@inputs = data.table(entities[, data.table::key(entities), with = FALSE])
        .Object@runinfo = data.table(entities[, data.table::key(entities), with = FALSE])        
        .Object@runinfo$outdir = paste(rootdir, task@name, entities[[data.table::key(entities)]], sep = '/') 
        .Object@runinfo$outdir = paste(rootdir, gsub('\\W', '.', task@name), entities[[data.table::key(entities)]], sep = '/') ## JEREMIAH
        .Object@runinfo[, status := 'ready']
        .Object@runinfo[, now := now]
        setkeyv(.Object@inputs, data.table::key(entities))
        setkeyv(.Object@outputs, data.table::key(entities))
        setkeyv(.Object@runinfo, data.table::key(entities))
        setkeyv(.Object@stamps, data.table::key(entities))
        
        if (length(ann.args)>0)
            {
                 ### of course we are noting timestamps <now>
                 ### the inputs may be modified between now and run time
                 ### in which case the task output may appear to be falsely outdated
                 ### when we check in the future ..
                 ### but this (false positive outdating) is less dangerous than false-negative outdating
### ie if we thought that the task was up to date but in fact inputs have changed.
                if (!mock)
                    cat('Noting time stamps of inputs\n')
                 for (this.arg in names(ann.args))
                     {
                         this.ann = ann.args[[this.arg]]@arg
                         nfiles = length(setdiff(entities[[this.ann]], NA))
                         
                         if (is.list(entities[[this.ann]]))
                         {
                             warning(sprintf('Annotation column "%s" is a list and not a vector.  Flow currently only supports entities tables with vector columns.  Will default to taking first element of each list item.  Please check the %s column to ensure that the data is correct.', this.ann, this.ann))
                             entities[[this.ann]] = sapply(entities[[this.ann]], '[', 1)
                         }

                         .Object@inputs[[this.arg]] = entities[[this.ann]]
        #                 .Object@inputs[, eval(this.arg) := entities[[this.ann]]]
                         
                         if (.Object@task@args[[this.arg]]@path)
                         {
                           if (!mock)
                           {
                             fn = .Object@inputs[[this.arg]]
                             fn[nchar(fn)==0] = NA ## NA out blank paths
                             nfiles = sum(!is.na(fn))
                             cat('\tfor', this.arg, sprintf('(%s files)', nfiles), '\n')
                           }
                           
                           if (is.numeric(entities[[this.ann]]))
                             stop(sprintf('Numeric annotation %s provided as path argument %s to task - check task configuration or entities table', this.ann, this.arg)) 
                           
                           .Object@stamps[, eval(this.arg) := as.character(file.info(ifelse(is.na(entities[[this.ann]]), '', entities[[this.ann]]))$mtime)]
                                        #                                 if (!mock)
                                        #                                     print(this.arg)
                           
                           cmd = paste(this.arg, ':= normalizePath(', this.arg, ')')
                           
                           .Object@inputs[!is.na(.Object@stamps[[this.arg]]), eval(parse(text = cmd))]
                           
                           if (!is.null(default(ann.args[[this.arg]])))
                             if (any(fenix <- !file.exists(ifelse(is.na(.Object@inputs[[this.arg]]), '', entities[[this.ann]]))))                                    
                                        #                                     if (any(fenix <- !file.exists(.Object@inputs[[this.arg]])))
                             {
                               .Object@inputs[fenix, eval(this.arg) := default(ann.args[[this.arg]])]
                               .Object@stamps[fenix, eval(this.arg) := as.character(Sys.time())]
                             }
                         }
                         else
                         {

                             nix = (nchar(str_trim(.Object@inputs[[this.arg]]))==0)
                             nix = nix[!is.na(nix)]
                             if (length(nix)>0)
                                 if (any(nix))
                                     .Object@inputs[nix, eval(this.arg) := NA]
                             
                             if (!is.null(default(ann.args[[this.arg]])))
                                 if (any(fenix <- is.na(.Object@inputs[[this.arg]])))
                                     .Object@inputs[[this.arg]][fenix] = default(ann.args[[this.arg]])
                             
                             .Object@stamps[, eval(this.arg) := ifelse(!is.na(.Object@inputs[[this.arg]]), as.character(Sys.time()), NA)]
                         }
                     }                                                     
             }

         if (length(lit.args)>0)
         {
                 .Object@inputs = cbind(.Object@inputs, as.data.table(as.data.frame(lit.args)[rep(1, nrow(entities)), , drop = FALSE]))

                 for (this.arg in names(lit.args))
                     {
                         if (.Object@task@args[[this.arg]]@path)
                             .Object@stamps[,
                                            eval(this.arg) :=
                                                as.character(file.info(lit.args[[this.arg]])$mtime)]
                         else
                             .Object@stamps[, eval(this.arg) := as.character(Sys.time())]
                     }                
             }

        if (ncol(.Object@stamps)>1) ## this should basically always be true
            {
                cols = setdiff(names(.Object@stamps), data.table::key(.Object))
                ready.mat = is.na(as.matrix(.Object@stamps[, cols, with = FALSE]))
                .Object@runinfo$status = ifelse(rowSums(ready.mat)>0, 'not ready', 'ready')                
                .Object@runinfo$status.info = ''
                if (any(ix <- .Object@runinfo$status != 'ready'))                    
                    .Object@runinfo$status.info[ix] = apply(ready.mat[ix, , drop = FALSE], 1,
                                                   function(x) paste(paste(colnames(ready.mat)[x], collapse = ','), 'not ready'))
                
                if (any(ix <- nchar(.Object@runinfo$status.info)>1))
                    warning(sprintf('missing annotations resulting causing %s jobs to be not ready.\n Breakdown of detailed statuses (with # entities with each specific status):\n\t%s',
                                  sum(status(.Object)!='ready'),
                                  paste(tabstring(table(.Object@runinfo$status.info[ix]), sep = ' '), collapse = ',')))
                
                
            }
        else
            {
                .Object@runinfo$status = 'ready'
                .Object@runinfo$status.info = ''
            }

        if (!mock)
            {
                cat('making output directories under', paste(.Object@rootdir, gsub('\\/', '.', task@name), sep = '/'), '\n')
                sapply(paste('mkdir -pv', .Object@runinfo$outdir), system)
#                system(paste('mkdir -p', paste(.Object@runinfo$outdir, collapse = ' ')))
            }
        
        ## instantiate commands for each row in entities table         
        .Object@runinfo[, cmd.og := sapply(1:nrow(entities), function(this.entity)
            {
                this.cmd = str_replace_all(module@cmd, stringr::fixed('<libdir>'), task@libdir)                
                for (k in 1:length(task@args))
                {
                        this.arg = names(task@args)[[k]]
                        this.val = .Object@inputs[this.entity, ][[gsub('\\-', '.', names(task@args)[k])]]
                        if (is.na(this.val))
                            this.val = ''
                        this.cmd = str_replace_all(this.cmd, stringr::fixed(paste('<', this.arg, '>', sep = '')), this.val)
                    }
                ## final cmd with all placeholders replaced
                return(this.cmd)
            })]

        if (is.null(mem))
            mem = NA
            
        if (!is.null(mem))
            mem = cbind(1:nrow(entities), mem)[,2]

        setkeyv(.Object@stamps, data.table::key(entities))
        
        .Object@runinfo[, queue := queue]
        .Object@runinfo[, mem := ifelse(is.na(mem), task@mem, mem)]
        .Object@runinfo[, nice := nice]
        .Object@runinfo[, cores := cores]

        ## append time call and input / output redirects to local function calls
        .Object@runinfo[, stderr := paste(outdir, '/', task@name, '.', entities[[data.table::key(entities)]], '.bsub.err', sep = '')]
        .Object@runinfo[, stdout := paste(outdir, '/',  task@name, '.', entities[[data.table::key(entities)]], '.bsub.out', sep = '')]

        .Object@runinfo[, jname := .jname(outdir, task@name, entities[[data.table::key(entities)]])]

#        if (!mock)
        .Object@runinfo = .update_cmd(.Object) ## updates BCMD / CMD


        ## populate output collection
        if (length(task@outputs)>0)
            if (sum(sapply(task@outputs, is, 'FlowOutput')>0))
                {
                    if (!mock)
                        cat('initializing output annotations\n')
                    for (nm in sapply(task@outputs, function(x) x@name))
                        .Object@outputs[, eval(nm) := as.character(NA)]
                }
        
        ## dump out rds of Job file for each entity, so we can reconstruct jobs, check outputs, compare
        ## time stamps etc later on
        if (!mock)
            {
                cat('Dumping out', length(.Object), 'Job.rds files to subdirectories of', .Object@rootdir, '\n')
                sapply(1:length(.Object), function(x) saveRDS(.Object[x], paste(.Object@runinfo$outdir[x], 'Job.rds', sep = '/')))

                cache(.Object)
                ids = ids(.Object)
                update(.Object, check.inputs = TRUE)
                path = paste(.Object@rootdir, '/', task(.Object)@name, '.rds', sep = '')
#                cat('Refreshing object from', path, '\n')
                .Object = readRDS(path)[ids, id = TRUE]
            }
        return(.Object)
    })

.jname = function(outdir, name, ids) paste(outdir, '/', name, '.', ids, sep = '')

.update_cmd = function(.Object)
    {
        ## utility func for instantiation of Job and modifying memory
        .cmd2bcmd = function(cmd, outdir, name, ids, queue, mem, cores) bsub_cmd(paste('touch ', outdir, '/started; ', cmd, ';', sep = ''), queue = queue, mem = mem, mc.cores = cores, cwd = outdir, jname = .jname(outdir, name, ids), jlabel = .jname(outdir, name, ids))
        .cmd2qcmd = function(cmd, outdir, name, ids, queue, mem, cores, now) qsub_cmd(cmd, queue = queue, mem = mem, mc.cores = cores, cwd = outdir, jname = paste('job', name, ids, sep = '.'), jlabel = paste('job', name, ids, sep = '.'), now = now)
        
        .Object@runinfo[, bcmd := '']
        ix = which(status(.Object) != 'not ready')
        .Object@runinfo[ix, bcmd := .cmd2bcmd(cmd.og, outdir, .Object@task@name, ids(.Object)[ix], queue, mem, cores)]
        .Object@runinfo[, cmd := '']
        .Object@runinfo[, cmd.quiet := '']
        .Object@runinfo[ix, cmd := paste('flow_go=$( pwd ); cd ', outdir, ';touch ', outdir, '/started; ', ifelse(nice, '(ionice -c2 -n7 nice ', ''), '/usr/bin/time -v ', cmd.og, ' ) 2>&1 | tee ', stdout, '; cp ', stdout, ' ', stderr, ';cd $flow_go',  sep = '')]
        .Object@runinfo[ix, cmd.quiet := paste('flow_go=$( pwd ); cd ', outdir, ';touch ', outdir, '/started; ', ifelse(nice, 'ionice -c2 -n7 nice ', ''), '/usr/bin/time -v ', cmd.og, ' &> ', stdout, '; cp ', stdout, ' ', stderr, ';cd $flow_go',  sep = '')]

        .Object@runinfo$cmd.path = paste(outdir(.Object), '/', names(outdir(.Object)), '.cmd.sh', sep = '')

        ## write cmd.og to file for qsub command
#        .Object@runinfo[, mapply(function(text, path) writeLines(text, path), paste('flow_go=$( pwd ); cd ', outdir, ';touch ', outdir, '/started; /usr/bin/time -v ', cmd.og, ' &> ', stdout, '; cp ', stdout, ' ', stderr, ';cd $flow_go',  sep = ''), cmd.path)] ## writes cmd to path        
                                        #        .Object@runinfo[, mapply(function(text, path) writeLines(text, path), paste('echo "FLOW.SGE.JOBID=$JOB_ID"; cd ', outdir, ';touch ', outdir, '/started; ~/Software/time/time -v ', cmd.og, '; cp ', stdout, ' ', stderr, sep = ''), cmd.path)] ## writes cmd to path

        .Object@runinfo[, mapply(function(text, path) writeLines(text, path), cmd, cmd.path)] ## writes cmd to path        
        .Object@runinfo[ix, qcmd := .cmd2qcmd(cmd.path, outdir, .Object@task@name, ids(.Object)[ix], queue, mem, cores, now = now)]
        
        return(.Object@runinfo)
    }


#' @name Job
#' @title Constructs a job object from a Task and keyed data.table of one or more entities. 
#'
#' Job instantiation combines the Task configuration with job specific info to create $cmd, $bcmd, $qcmd to run job locally and submit to LSF / SGE and output
#' to task / entity specific output directories,  The job object can be polled to examine job status, and run or re-run specific jobs. 
#'
#' A Job object is instantiated from the combination of a text .task file containing the task configuration
#' and a keyed data.table of entities (eg samples, pairs, individuals)
#' Instantation involves populating the task with all the relevant columns of the entitie table,
#' 
#' 
#' Job instantiation requires a Task object or path to a .task file as input + and keyed "entities" data table with all the necessary input columns
#' required by the task configuration.
#' If any of these columns don't exist then the Job object will fail to be instantiated.  Optional input rootdir specifies root directory of all job specific  output directories. 
#' All outputs for jobs will be written to directories that are named by the task name / row key.   These directories are created
#' at time of object instantiation.  An .rds of the job object is stored in a standard location in the rootdir (with name TASKNAME.rds)
#'
#' @param task task config (.task file) or Task object
#' @param entities keyed data table of entities that contain annotations which task will be drawing from
#' @param rootdir the root directory under which Task specific output will be placed (default ./Flow)
#' @param queue which queue to direct jobs to (should be compatible with local LSF or SGE HPC)
#' @param mem memory limit to put on jobs (in GB)
#' @param nice whether to nice jobs when running locally (default = TRUE)
#' @param cores how many cores to run jobs with
#' @param mock boolean, if FALSE will not create subdirectories
#' @author Marcin Imielinski
#' @export
Job = function(
    task, ##Task wrapping around an Module expecting literal and annotation arguments
    entities, ## keyed data.table, key will determine id of outgoing jobs, columns of table used to populate task
    rootdir = './Flow/',
    queue = as.character(NA),
    mem = NULL,
    nice = NULL,
    cores = 1,
    mock = FALSE, ...) {
    new('Job', task = task, entities = entities, rootdir = rootdir,
        queue = queue, nice = nice, mem = mem, cores = cores, mock = mock, ...)
}


#' @name c
#' @title Concatenates Job objects built from the same Task
#' 
#'
#' @exportMethod c
#' @export
#' @author Marcin Imielinski
setMethod('c', 'Job', function(x, ...)
    {
        if (!.hasSlot(x, 'entities'))
            stop('older version of Flow object does not support concatenation')
        
        obj = c(list(x), list(...))
        same.same = all(sapply(obj, function(x) identical(x@task, obj[[1]]@task)))
        if (!all(same.same))
            stop('trying to conatenate incompatible Job objects: can only concatenate Jobs built from same Task')

        ukey = unique(sapply(obj, function(x) data.table::key(x@inputs)))
        if (length(ukey)>1)
            stop('trying to concatenate incompatible Job objects built from different keys')
        
        urootdir = unique(sapply(obj, function(x) x@rootdir))
        if (length(urootdir)>1)
            warning('concatenating  Job objects with different rootdirs, choosing first for output')
        
        icol = names(obj[[1]]@inputs)
        ocol = names(obj[[1]]@outputs)
        scol = names(obj[[1]]@stamps)
        rcol = names(obj[[1]]@runinfo)
        ecol = names(obj[[1]]@entities)
                
        inputs = rbindlist(lapply(obj, function(x) x@inputs[, icol, with = FALSE]))
        outputs = rbindlist(lapply(obj, function(x) x@outputs[, ocol, with = FALSE]))
        runinfo = rbindlist(lapply(obj, function(x) x@runinfo[, rcol, with = FALSE]))
        stamps = rbindlist(lapply(obj, function(x) x@stamps[, scol, with = FALSE]))
        entities = rbindlist(lapply(obj, function(x) x@entities[, ecol, with = FALSE]))

        setkeyv(entities, ukey)
        setkeyv(inputs, ukey)
        setkeyv(outputs, ukey)
        setkeyv(runinfo, ukey)
        setkeyv(stamps, ukey)
        
        return(Job(obj[[1]]@task, entities = entities, rootdir = urootdir[1], queue = runinfo$queue, mem = runinfo$mem, nice = runinfo$nice, cores = runinfo$cores, now = runinfo$now, mock = TRUE)) 
    })

#' @name purge
#' @exportMethod purge
#' @export
setGeneric('purge', function(object, ...) {standardGeneric('purge')})


#' @name purge
#' @title purges all output / run directories associated with this Job object (be careful)
#' 
#' @exportMethod purge
#' @export
#' @author Marcin Imielinski
setMethod('purge', 'Job', function(object, check.inputs = TRUE)
    {
        cat('About to remove all files and directories associated with this Job object in', paste(object@rootdir, task(object)@name, sep = '/'), '\nGiving you a moment to think ... ')
        Sys.sleep(5)
        cat('OK here we go .. \n')
        sapply(outdir(object), function(x) system(paste('rm -r', x)))
        cat('Regenerating fresh output directories\n')
        sapply(outdir(object), function(x) system(paste('mkdir -p', x)))
        sapply(1:length(object), function(x) saveRDS(object[x], paste(outdir(object)[x], 'Job.rds', sep = '/')))
        cat('Done\n')
    })


#' @name Job-class
#' @rdname Job-class
#' @exportMethod update
#' @export
setGeneric('update', function(object, ...) {standardGeneric('update')})


#' @name update
#' @title updates Job status
#'
#' This method checks inputs for updates to input files, polls output directories for the presence output files, checks stdout output logs, among
#' other things to assess job status.  Once finished, it caches results to the .rds object associated with this Job. 
#' 
#' @exportMethod update
#' @export
#' @author Marcin Imielinski
setMethod('update', 'Job', function(object, check.inputs = TRUE, mc.cores = 1, cache.object = TRUE, print.status = TRUE)
{
        ## for every output, apply regexp in files of outdir to search for files
        status.info = rep('', length(object))
        status = rep('ready', length(object))

        new.object = object
        ids = new.object@outputs[[data.table::key(new.object@outputs)]]
        
        st = file.info(paste(outdir(new.object), 'started', sep = '/'))
        en = file.info(paste(outdir(new.object), 'failed', sep = '/'))
        rep = report(new.object, mc.cores = mc.cores)
        status = ifelse(!is.na(st$mtime),
            ifelse(!is.na(rep$success),
                   ifelse(rep$success, 'completed', 'failed'),
                   ifelse(!is.na(st$mtime), 'running', 'ready')), 'ready')

        
        if (length(new.object@task@outputs)>0) ## check output args if they exist       
            if (sum(sapply(new.object@task@outputs, is, 'FlowOutput')))
                {
                    outkeys = sapply(new.object@task@outputs, function(x) x@name)
                    outre = sapply(new.object@task@outputs, function(x) x@pattern)
                    for (id in ids)
                        {
                            files = dir(new.object@runinfo[id, outdir])
                            names(files) = paste(new.object@runinfo[id, outdir], files, sep = '/')                           
                            for (k in 1:length(outkeys))
                                new.object@outputs[id, eval(outkeys[k]) := names(files)[grep(outre[k], files)][1]]
                        }

                    out.status = !is.na(new.object@outputs[, outkeys, with = FALSE])
                    has.out = rowSums(out.status)>0
                    missing.out = rowSums(!out.status)>0
                    status = ifelse(status=='completed',
                        ifelse(missing.out, 'completed; some outputs missing', 'completed'),
                        ifelse(has.out, paste(status, 'and some outputs present'), status))
                }
        
        ## determine ready / not ready / outdated status based on the existence of
        ## of file names
        args = new.object@task@args

        if (length(args)>0)
            if (check.inputs)
            {
              output.date = as.POSIXct(file.info(new.object@runinfo$stdout)$mtime)
              
              outdated = matrix(FALSE, nrow = length(new.object), ncol = length(args), dimnames = list(ids, names(args)))
              cat('Checking input date stamps\n')
              for (this.arg in names(args))
              {
                if (args[[this.arg]]@path)
                {
                  fn = new.object@inputs[[this.arg]]
                  fn[nchar(fn)==0] = NA ## NA out blank paths
                  nfiles = sum(!is.na(fn))
                  cat('\tfor', this.arg, sprintf('(%s files)', nfiles), '\n')
                  fe = file.exists(fn)
                  old.date = as.POSIXct(new.object@stamps[[this.arg]])
                  if (any(fe))
                    {
                      if (any(ix<-is.na(old.date))) ## if for some reason blank, set to some time in the far future
                        old.date[ix] = Sys.time()+5e9
                      if (is(args[[this.arg]], 'FlowLiteral') & args[[this.arg]]@path)
                        outdated[, this.arg] =
                                        #                                          as.POSIXct(as.character(file.info(args[[this.arg]]@arg)$mtime))>old.date &
                          ifelse(is.na(output.date), FALSE, as.POSIXct(as.character(file.info(new.object@inputs[[this.arg]])$mtime))>output.date)
                      else if (is(args[[this.arg]], 'FlowAnnotation') & args[[this.arg]]@path)
                        outdated[, this.arg] =
                                        #                                         as.POSIXct(as.character(file.info(new.object@inputs[[this.arg]])$mtime))>old.date &
                          ifelse(is.na(output.date), FALSE, as.POSIXct(as.character(file.info(new.object@inputs[[this.arg]])$mtime))>output.date)
                      else
                        outdated[, this.arg] = FALSE
                    }

                  if (any(!fe))
                    outdated[!fe, this.arg] = NA
                }
              }


              status = ifelse(rowSums(outdated, na.rm = TRUE)>0, 'outdated', status)
              status.info = paste(status.info, apply(outdated, 1,
                                                     function(x) if (length(which(x))>0) paste('Updates in', paste(colnames(outdated)[which(x)], collapse = ', '))
                                                                 else ''))

                notready = rowSums(is.na(outdated))>0
                if (any(notready))
                {
                    status[notready] = 'not ready'
                    status.info[notready] = paste(status.info[notready], apply(is.na(outdated[notready, , drop = FALSE]), 1,
                                           function(x) paste(paste(colnames(outdated)[which(x)], collapse = ', '), 'not ready')))
                }
            }
    
        new.object@runinfo$status = status
        new.object@runinfo$status.info = str_trim(status.info)
        new.object@runinfo = .update_cmd(new.object)

        if (cache.object)
            cache(new.object)

        if (print.status)
            print(table(status(new.object)))
        ## weird R voodoo to modify object in place
        eval(
              eval(
                       substitute(
                                   expression(object <<- new.object)
                                ,env=parent.frame(1) )
                    )
        )
        cat('')
    })


setGeneric('default', function(.Object) standardGeneric('default'))

setMethod('default', 'FlowAnnotation', function(.Object)
    {
        if (.hasSlot(.Object, 'default'))
            if (length(.Object@default)==1)
                .Object@default
            else
                NULL
        else
           NULL
    })

setMethod('default', 'FlowLiteral', function(.Object)
    {
        NULL
    })



#' @export
setGeneric('cache', function(object, ...) {standardGeneric('cache')})


#' @name cache
#' @title caches .rds copy of job object to standard location (TASK.NAME.rds in Job output root directory)
#' 
#' @exportMethod cache
#' @export
#' @author Marcin Imielinski
setMethod('cache', 'Job', function(object, verbose = TRUE)
    {
        path = paste(object@rootdir, '/', task(object)@name, '.rds', sep = '')

        if (file.exists(path))
            {
                old.cached = readRDS(path)                
                if (length(others <- setdiff(ids(old.cached), ids(object)))>0)
                    {
                        new.ids = union(ids(old.cached), ids(object))
                        object = tryCatch(c(old.cached[others, id = TRUE], object), error = function(e) NULL)
                        if (is.null(object))
                            stop(sprintf('Error merging object with existing object in %s, please check format of stored object and fix or delete the corresponding .rds file', path))
                        object = object[new.ids, id = TRUE]
                cat('Caching object to', path, 'after merging with', length(others), 'additional entities \n')
                    }
                else
                    cat('Caching object to', path, '\n')
            }
        else
            cat('Caching object to', path, '\n')
            
        saveRDS(object, path)
    })


#' @export
setGeneric('chmod', function(object, ...) {standardGeneric('chmod')})


#' @name chmod
#' @title chmods Job output directories
#'
#' Chmods outputs of jobs (by default to 775)
#' 
#' @exportMethod chmod
#' @export
#' @author Marcin Imielinski
setMethod('chmod', 'Job', function(object, mode = '775', verbose = FALSE, mc.cores = 1)
{
  if (verbose)
    verbose = '-v'
  else
    verbose = ''

  mclapply(outdir(object), function(x) system(paste('chmod -R ', verbose, mode, x)),
                                              mc.cores = mc.cores)
})


#' @export
setGeneric('refresh', function(object, ...) {standardGeneric('refresh')})


#' @name refresh
#' @title refreshes .rds copy of job object from standard location (TASK.NAME.rds in Job output root directory)
#' 
#' @exportMethod cache
#' @export
#' @author Marcin Imielinski
setMethod('refresh', 'Job', function(object, verbose = TRUE)
    {
        path = paste(object@rootdir, '/', task(object)@name, '.rds', sep = '')
        cat('Refreshing object from', path, '\n')
        new.object = readRDS(path)
        eval(
            eval(
                substitute(
                    expression(object <<- new.object)
                   ,env=parent.frame(1) )
                )
            )
        cat('')
    })

#' @name [
#' @title caches .rds copy of job object to standard location (TASK.NAME.rds in Job output root directory)
#' @description
#' 
#' Subsetting Job object, can use index or character.  If character will act as a grep
#' of status, or status.info if status.info = TRUE
#' 
#' @exportMethod [
#' @export
#' @author Marcin Imielinski
setMethod('[', 'Job', function(x, i, id = FALSE)
    {
        if (!id & is.character(i))
            i = grep(i, status(x), ignore.case = TRUE)        
        else
            if (is.logical(i))
                i = which(i)

        if (length(i) ==0)
          i = 0 ## data.table-ese for empty data.table, otherwise NULL

        id = key(x)
        x@runinfo = x@runinfo[i, ]
        x@inputs = x@inputs[i, ]
        x@outputs = x@outputs[i, ]
        x@stamps = x@stamps[i, ]
        if (.hasSlot(x, 'entities'))
        {
            setkeyv(x@entities, id)
            x@entities = x@entities[i, ]
            setkeyv(x@entities, id)
        }
        
        ## some kind of data.table bug where key gets lost every w subsetting once in a while ... #CHECK
        setkeyv(x@runinfo, id)
        setkeyv(x@inputs, id)
        setkeyv(x@outputs, id)
        setkeyv(x@stamps, id)
        return(x)
    })


#' Replacing a subset of Job object with another (compatible Job object) 
#'
#' @name [<-
#' @aliases [,Job-class
#' @rdname Job-class
#' @docType methods
#' @exportMethod [<-
#' @export

#' @name [<-
#' @title Replaces one or more Job items
#' @description
#' 
#' Replacing a subset of Job object with other compatibler Job items.  "compatible" job items
#' will be from the identical task and will contain entities already contained in the Job
#' object. 
#' 
#' @exportMethod [<-
#' @export
#' @author Marcin Imielinski
setMethod('[<-', 'Job', function(x, i, value, id = FALSE)
    {          
        if (!is(value, 'Job'))
            stop('Must replace with Job object')

        if (!(identical(x@task, value@task)))
            stop('Can only sub in Job having identical task')

        if (!(identical(names(x@inputs), names(value@inputs))))
            stop('Can only sub in Job having identical task')

        if (!id & is.character(i))
            i = grep(i, status(x), ignore.case = TRUE)        
        else
            if (is.logical(i))
                i = which(i)

        if (ncol(value@runinfo) != ncol(x@runinfo))
            stop('Either object or replacement value are corrupt: please try regenerating Job object')
        x@runinfo[i,] = value@runinfo
        x@inputs[i,] = value@inputs
        x@outputs[i, ] = value@outputs
        x@stamps[i, ] = value@stamps
        if (.hasSlot(x, 'entities') & .hasSlot(value, 'entities'))
            x@entities[i, ] = value@entities
        return(x)
    })



#' @name Job-class
#' @rdname Job-class
#' @exportMethod module
#' @export
setGeneric('module', function(.Object) { standardGeneric('module')})


#' @name module
#' @title Gets module associated with a Job
#' @exportMethod module
#' @export
#' @author Marcin Imielinski
setMethod('module', 'Job', function(.Object)
    {
        .Object@task@module
    })



#' @name module-task
#' @title Gets module associated with a Task
#' @description
#' getting module associated with Task object
#'
#' @exportMethod module
#' @export
#' @author Marcin Imielinski
setMethod('module', 'Task', function(.Object)
    {
        .Object@module
    })


#' @name Job-class
#' @rdname Job-class
#' @exportMethod task
#' @export
setGeneric('task', function(.Object) {standardGeneric('task')})


#' @name task
#' @title gets Task associated with Job object
#' 
#' @exportMethod task
#' @author Marcin Imielinski
#' @export
setMethod('task', 'Job', function(.Object)
    {
        .Object@task
    })

#' @name Job-class
#' @rdname Job-class
#' @exportMethod mem
#' @export
setGeneric('mem', function(.Object) {standardGeneric('mem')})

#' @name mem
#' @title Gets max memory associated with the jobs in the Job object
#' @description
#' getting max mem (in GB) associated with Job object
#'
#' @exportMethod mem
#' @export
#' @author Marcin Imielinski
setMethod('mem', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, mem], names = .Object@runinfo[[key(.Object@runinfo)]])
    })



#' @export
setGeneric('mem<-', function(.Object, value) {standardGeneric('mem<-')})



#' @name mem<-
#' @title Sets max memory associated with the jobs in the Job object
#' @description
#' Setting max mem (in GB) associated with Job object
#'
#' @exportMethod mem<-
#' @export
#' @author Marcin Imielinski
setReplaceMethod('mem', 'Job', function(.Object, value)
                 {
                     .Object@runinfo[, mem := value]                     
                     .Object@runinfo = .update_cmd(.Object)
                     return(.Object)
                 })

#' @name Job-class
#' @rdname Job-class
#' @exportMethod queue
#' @export
setGeneric('queue', function(.Object) {standardGeneric('queue')})

#' @name queue
#' @title Gets queue associated with the jobs in the Job object
#' @description
#' Getting LSF / SGE queue associated with Job object
#'
#' @exportMethod queue
#' @export
#' @author Marcin Imielinski
setMethod('queue', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, queue], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })


#' @export
setGeneric('queue<-', function(.Object, value) {standardGeneric('queue<-')})

#' @name queue<-
#' @title Sets queue associated with the jobs in the Job object
#' @description
#' Setting LSF / SGE queue associated with Job object
#'
#' @exportMethod queue<-
#' @export
#' @author Marcin Imielinski
setReplaceMethod('queue', 'Job', function(.Object, value)
                 {
                     .Object@runinfo[, queue := value]
                     .Object@runinfo = .update_cmd(.Object)
                     return(.Object)
                 })

######
#' @name Job-class
#' @rdname Job-class
#' @exportMethod cores
#' @export
setGeneric('cores', function(.Object) {standardGeneric('cores')})

#' @name cores
#' @title Gets cores associated with the jobs in the Job object
#' @description
#' Getting LSF / SGE cores associated with Job object
#'
#' @exportMethod cores
#' @export
#' @author Marcin Imielinski
setMethod('cores', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, cores], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })


#' @export
setGeneric('cores<-', function(.Object, value) {standardGeneric('cores<-')})

#' @name cores<-
#' @title Sets cores associated with the jobs in the Job object
#' @description
#' Setting LSF / SGE cores associated with Job object
#'
#' @exportMethod cores<-
#' @export
#' @author Marcin Imielinski
setReplaceMethod('cores', 'Job', function(.Object, value)
                 {
                     .Object@runinfo[, cores := value]
                     .Object@runinfo = .update_cmd(.Object)
                     return(.Object)
                 })




#' @name Job-class
#' @rdname Job-class
#' @exportMethod now
#' @export
setGeneric('now', function(.Object) {standardGeneric('now')})


#' @name now
#' @title Toggles whether to force jobs to run now
#' @description
#'
#' @exportMethod now
#' @export
#' @author Marcin Imielinski
setMethod('now', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, now], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })





#' @export
setGeneric('now<-', function(.Object, value) {standardGeneric('now<-')})

#' @name now<-
#' @title Sets now param associated with the jobs in the Job object
#' @description
#' Setting LSF / SGE now associated with Job object
#'
#' @exportMethod now<-
#' @export
#' @author Marcin Imielinski
setReplaceMethod('now', 'Job', function(.Object, value)
                 {
                     .Object@runinfo[, now := value]
                     .Object@runinfo = .update_cmd(.Object)
                     return(.Object)
                 })


#' @export
setGeneric('cmd', function(.Object, ...) {standardGeneric('cmd')})


#' @name cmd
#' @title Gets vector of command line calls associated with the jobs in the Job object
#' @exportMethod cmd
#' @export
#' @author Marcin Imielinski
setMethod('cmd', 'Job', function(.Object, all = FALSE, quiet = TRUE)
    {
        if (!all)
            {
                ix = runinfo(.Object)[, which(!(status %in% c('completed', 'not ready')))]
                if (quiet)
                    if (is.null(.Object@runinfo$cmd.quiet))
                        structure(.Object@runinfo[, cmd[ix]], names = ids(.Object)[ix])
                    else
                        structure(.Object@runinfo[, cmd.quiet[ix]], names = ids(.Object)[ix])
                else
                    structure(.Object@runinfo[, cmd[ix]], names = ids(.Object)[ix])
                    
            }
        else
            {
                if (quiet)
                    structure(.Object@runinfo[, cmd.quiet], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
                else
                    structure(.Object@runinfo[, cmd], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
            }
    })


#' @export
setGeneric('run', function(.Object, ...) {standardGeneric('run')})


#' @name run
#' @title Runs jobs as local commands from shell via R system calls
#' @exportMethod run
#' @export
#' @param mc.cores Number of parallel cores to run jobs with (=1)
#' @author Marcin Imielinski
setMethod('run', 'Job', function(.Object, mc.cores = 1, all = FALSE, quiet = TRUE)
    {
        require(parallel)


        cmds = cmd(.Object, quiet = quiet, all = all)
        if (length(cmds)==0)
            {
                cat('No jobs to run\n')
                return()                
            }
        if (is.null(names(cmds)))
            names(cmds) = ids(.Object)        
                                        #        mclapply(names(cmd(.Object, quiet = quiet)), function(x)
        mclapply(names(cmds), function(x)
            {
                cat('Starting', task(.Object)@name, 'on entity',  x, '\n')
                                        #                system(cmd(.Object, quiet = quiet)[x])
                system(cmds[x])
            }, mc.cores = mc.cores)
    })

#' @export
setGeneric('bcmd', function(.Object, ...) {standardGeneric('bcmd')})

#' @name bcmd
#' @title Returns vector of bsub (ie LSF) commands associated with this Job object
#' @exportMethod bcmd
#' @export
#' @param all logical flag whether to run all jobs (including completed)
#' @author Marcin Imielinski
setMethod('bcmd', 'Job', function(.Object, all = FALSE)
    {
        if (!all)
            {
                ix = runinfo(.Object)[, which(!(status %in% c('completed', 'not ready')))]
                structure(.Object@runinfo[, bcmd[ix]], names = ids(.Object)[ix])
            }
        else
            structure(.Object@runinfo[, bcmd], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })

#' @export
setGeneric('brun', function(.Object, ...) {standardGeneric('brun')})


#' @name brun
#' @title Runs jobs on LSF
#' @exportMethod brun
#' @export
#' @param mc.cores Number of parallel cores to run jobs with (=1)
#' @author Marcin Imielinski
setMethod('brun', 'Job', function(.Object, mc.cores = 1, all = FALSE)
    {
        mclapply(bcmd(.Object, all = all), system, mc.cores = mc.cores)
    })


#' @export
setGeneric('qcmd', function(.Object, ...) {standardGeneric('qcmd')})

#' @name qcmd
#' @title Returns vector of qsub (ie SGE) commands associated with this Job object
#' @exportMethod qcmd
#' @export
#' @param all logical flag whether to run all jobs (including completed)
#' @author Marcin Imielinski
setMethod('qcmd', 'Job', function(.Object, all = FALSE)
{
    if (!('qcmd' %in% names(.Object@runinfo)))
    {
        warning('Job object corrupted: qcmd missing')
        .Object@runinfo = .update_cmd(.Object)
    }
        
        if (!all)
        {
            ix = runinfo(.Object)[, which(!(status %in% c('completed', 'not ready')))]
            structure(.Object@runinfo[, qcmd[ix]], names = ids(.Object)[ix])
        }
        else
            structure(.Object@runinfo[, qcmd], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })

#' @export
setGeneric('qrun', function(.Object, ...) {standardGeneric('qrun')})
#' @name qrun
#' @title Runs jobs on LSF
#' @exportMethod qrun
#' @export
#' @param mc.cores Number of parallel cores to run jobs with (=1)
#' @author Marcin Imielinski
setMethod('qrun', 'Job', function(.Object, mc.cores = 1, all = FALSE)
    {
        qcmds = qcmd(.Object, all = all)
        res = lapply(qcmds, function(x) {p = pipe(x); out = readLines(p); close(p); return(out)})
        jobids = sapply(res, function(x) gsub('Your job (\\d+) .*', '\\1', x))
        mapply(function(d,j) writeLines(j, paste0(d,'/sge.jobid')), outdir(.Object)[names(qcmds)], jobids) ## save last jobids
        writeLines(paste('Deploying', jobids, 'for entity', ids(.Object)))
    })



#' @export
setGeneric('outdir', function(.Object) {standardGeneric('outdir')})

#' @name outdir
#' @title Gets vector of output directories associated with the jobs in this Job object
#' @exportMethod outdir
#' @export
#' @author Marcin Imielinski
setMethod('outdir', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, outdir], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })


#' @export
setGeneric('dirs', function(.Object, ...) {standardGeneric('dirs')})


#' @name dirs
#' @title returns list of character vectors containing the contents of directory associated with these jobs
#' @exportMethod dirs
#' @export
#' @author Marcin Imielinski
setMethod('dirs', 'Job', function(.Object, pattern = NULL, full = TRUE, ...)
    {
        out = lapply(outdir(.Object), dir, pattern = pattern, full = full, ...)
        return(out)
    })

#' @export
setGeneric('err', function(.Object) {standardGeneric('err')})

#' @name err
#' @title returns vector of paths to stderr files associated with this job object
#' @exportMethod err
#' @export
#' @author Marcin Imielinski
setMethod('err', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, stderr], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })


#' @export
setGeneric('out', function(.Object) {standardGeneric('out')})

#' @name out
#' @title returns vector of paths to stdout files associated with this job object
#' @exportMethod out
#' @export
#' @author Marcin Imielinski
setMethod('out', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, stdout], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])
    })


#' @export
setGeneric('outputs', function(.Object) {standardGeneric('outputs')})


#' @name outputs
#' @title returns keyed data.table of output entity annotations generated from the jobs in this Job object.  Annotations will be NA for non-completed jobs.
#' @exportMethod outputs
#' @export
#' @author Marcin Imielinski
setMethod('outputs', 'Job', function(.Object)
    {
        copy(.Object@outputs)
    })


#' @export
setGeneric('runinfo', function(.Object) {standardGeneric('runinfo')})

#' getting output data.table associated with Job object
#'
#' @name runinfo
#' @aliases [,Job-class
#' @rdname Job-class
#' @docType methods
#' @exportMethod runinfo
#' @export

#' @name runinfo
#' @title returns keyed data.table of runtime parameters associated with this this job object, which will include $cmd, $bcmd, $mem, and $queue, among others.
#' @exportMethod runinfo
#' @export
#' @author Marcin Imielinski
setMethod('runinfo', 'Job', function(.Object)
    {
        copy(.Object@runinfo)
    })

                    
#' @export
setGeneric('inputs', function(.Object) {standardGeneric('inputs')})

#' @name inputs
#' @title returns keyed data.table of inputs annotation associated with the jobs in this Job object
#' @export
#' @author Marcin Imielinski
setMethod('inputs', 'Job', function(.Object)
    {
        copy(.Object@inputs)
    })


#' @export
setGeneric('entities', function(.Object) {standardGeneric('entities')})

#' @name entities
#' @title returns original keyed data.table of entities associated with this Job object
#' @export
#' @author Marcin Imielinski
setMethod('entities', 'Job', function(.Object)
    {
        if (.hasSlot(.Object, 'entities'))
            copy(.Object@entities)
        else
            data.table(fsno[, 1, with = FALSE])

    })


#' @export
setGeneric('status', function(.Object) {standardGeneric('status')})


#' @name status
#' @title Returns character vector of status of jobs in this Job object.  Regexps on these statuses can be also used as character inputs to [] accessor
#' @exportMethod status
#' @export
#' @author Marcin Imielinski
setMethod('status', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, status], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])        
    })

                                                                                                                                                                                                      
#' @export
setGeneric('jname', function(.Object) {standardGeneric('jname')})


#' @name jname
#' @title Returns character vector of job names associated with this Job object.  
#' @exportMethod jname
#' @export
#' @author Marcin Imielinski
setMethod('jname', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, jname], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])        
    })

#' @export
setGeneric('bjobs', function(.Object, ...) {standardGeneric('bjobs')})


#' @name bjobs
#' @title Tracks down any LSF jobs associated with this Job object using a bjobs query (warning can be slow)
#' @exportMethod jname
#' @export
#' @author Marcin Imielinski
setMethod('bjobs', 'Job', function(.Object, mc.cores = 1)
    {
        res = lapply(paste('bjobs -J', jname(.Object), '2> /dev/null'),
            function(x) {p = pipe(x); out = readLines(p); close(p); return(out)})
        res = mclapply(res, function(x) if (length(x)>0) sapply(strsplit(x, '[ ]+', perl = TRUE), function(y) y[c(1,3)])[,2],
                       mc.cores = mc.cores)       
        names(res) = ids(.Object)
        out = data.table(runinfo(.Object)[, key(.Object), with = FALSE])
        out$JOBID = as.character(NA)
        out$STATUS = as.character(NA)
        if (any(ix <- sapply(res, length)>0))
            out[ix, ':='(JOBID = sapply(res[ix], function(x) x[1]), STATUS = sapply(res[ix], function(x) x[2]))]
        return(out)            
    })



#' @export
setGeneric('qjobs', function(.Object, ...) {standardGeneric('qjobs')})


#' @name qjobs
#' @title Tracks down any sGE jobs associated with this Job object using a bjobs query (warning can be slow)
#' @exportMethod jname
#' @export
#' @author Marcin Imielinski
setMethod('qjobs', 'Job', function(.Object)
    {
                fn.jids = sapply(outdir(.Object), function(x) paste(x, 'sge.jobid', sep = '/'))
        ix = file.exists(fn.jids)
        out1 = out2 = NULL
        nms = c('jobid','prior','ntckt','name','user','project','department','state','cpu','mem','io','tckts','ovrts','otckt','ftckt','stckt','share','queue','slots')
#        nms = c('jobid', 'prior', 'name', 'user', 'state', 'start.sumit.at', 'queue', 'slots', 'taskid')
        out = runinfo(.Object)[, key(.Object), with = FALSE]
        for (nm in nms)
            out[[nm]] = as.character(NA)
        if (any(ix))
            {
                jids = rep(NA, length(.Object))
                jids[ix] = sapply(fn.jids[ix], function(x) readLines(x)[1])
                p = pipe('qstat -ext')
                tab = strsplit(str_trim(readLines(p)), '\\s+')
                close(p)
                iix = sapply(tab, length)<=length(nms) & sapply(tab, length)>14
                if (length(tab)>0)
                    {
                        tab = lapply(tab, function(x) x[1:length(nms)])
                        tmp = as.data.table(matrix(unlist(tab[iix]), ncol = length(nms), byrow = TRUE))
                        setnames(tmp, nms)
                        setkey(tmp, jobid)
                        out = cbind(runinfo(.Object)[, key(.Object), with = FALSE], tmp[jids, ])
                        na = is.na(out$state)
                        if (any(na))
                            out$jobid[na] = NA
                    }
            }
                return(out)            
            })


#' @export
setGeneric('qkill', function(.Object, ...) {standardGeneric('qkill')})

#' @name qkill
#' @title Kills any running SGE jobs associated with this Job object
#' @exportMethod qkill
#' @export
#' @author Marcin Imielinski
setMethod('qkill', 'Job', function(.Object, jid = NULL)
    {
        
        qj = qjobs(.Object)
        ix = !is.na(qj$jobid)
        if (any(ix))
            system(paste(c('qdel', qj$jobid[ix]), collapse = ' '))
        else
            cat('No queued or running SGE jobs to kill\n')
    })


#' @export
setGeneric('bkill', function(.Object, ...) {standardGeneric('bkill')})

#' @name bkill
#' @title Kills any running LSF jobs associated with this Job object
#' @exportMethod bkill
#' @export
#' @author Marcin Imielinski
setMethod('bkill', 'Job', function(.Object, jid = NULL)
    {
        sapply(paste('bkill -J', jname(.Object), '2>/dev/null'), system)
        cat('')
    })

#' @export
setGeneric('status.info', function(.Object) {standardGeneric('status.info')})



#' @name status.info
#' @title returns vector of more detailed statuses associated with jobs in this Job object
#' @exportMethod status.info
#' @export
#' @author Marcin Imielinski
setMethod('status.info', 'Job', function(.Object)
    {
        structure(.Object@runinfo[, status.info], names = .Object@runinfo[[data.table::key(.Object@runinfo)]])        
    })

#' @export
setGeneric('ids', function(.Object) { standardGeneric('ids')})


#' @name ids
#' @title returns vector of entity ids associated with this Job objedct
#' @exportMethod ids
#' @export
#' @author Marcin Imielinskix
setMethod('ids', 'Job', function(.Object)
    {
        .Object@runinfo[[key(.Object)]]
    })

#' @name ids
#' @title returns vector of entity ids associated with this Job objedct
#' @exportMethod ids
#' @export
#' @author Marcin Imielinskix
setMethod('ids', 'Job', function(.Object)
    {
        .Object@runinfo[[key(.Object)]]
    })


#' @name length
#' @title returns length of this Job object
#' @exportMethod length
#' @export
#' @author Marcin Imielinski
setMethod('length', 'Job', function(x)
    {
              return(nrow(x@outputs))
    })
    

#' @export
setGeneric('key', function(x) {standardGeneric('key')})

#' @name key
#' @title Retrieves the entity key associated with this Job object
#' @exportMethod key
#' @export
#' @author Marcin Imielinski
setMethod('key', 'Job', function(x)
    {
        return(data.table::key(x@outputs))
    })
    
setMethod('show', 'Job', function(object)
    {
        .tabstring = function(tab, sep = ', ')
            return(paste(names(tab), ' (', tab, ')', sep = '', collapse = sep))

        if (length(object)==1)
            estring = ids(object)
        else
            estring = paste0(substr(paste(ids(object), collapse = ', '), 1, 20), '...')
        
        cat(sprintf('Job on %s entities (%s) with rootdir %s from task %s using module %s version %s\nJob status: %s\n', length(object), estring, object@rootdir, object@task@name, object@task@module@name, object@task@module@stamp, .tabstring(table(status(object)))))
    })



### 
### util functions
###
###

##################
# Makes bsub command that wraps shell command "cmd" to send to queue "queue"
# redirebmccting output / error etc streams to path prefixed by "jname",
# optional_args: maximum memory requirements "mem", "jlabel" job label
##################
bsub_cmd = function(cmd, queue, jname = NULL, jlabel=NULL, jgroup = NULL, mem=NULL, group = "cgafolk", cwd = NULL, mc.cores = NULL, deadline = F)
  {
    if (is.null(jname) & is.null(names(cmd)))
      jname = 'job'

    if (length(jname) != length(cmd))
      jname = rep(jname, length(cmd))
    
    if (!is.null(jname))
      names(cmd) = dedup(jname)    
                            
    qjname = paste( "\"", names(cmd), "\"", sep="" )
    qjout = paste( "\"", names(cmd), ".bsub.out", "\" ", sep="" )
    qjerr = paste( "\"", names(cmd), ".bsub.err", "\" ", sep="" )
    qjrout = paste( "\"", names(cmd), ".R.out", "\" ", sep="" )
    out_cmd = paste( "bsub -o ", qjout, " -e ",  qjerr, " -P ", group);
    out_cmd = paste(out_cmd, ifelse(is.na(queue), '', paste("-q ", queue)))
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -R \"rusage[mem=", mem, "]\" ", sep = "");
    if (!is.null(jlabel)) out_cmd = paste(out_cmd, " -J ", jlabel )
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) out_cmd = paste(out_cmd, " -cwd ", cwd )
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, sprintf(" -n %d,%d -R 'span[hosts=1]'", mc.cores, mc.cores))
    if (deadline) out_cmd = paste(out_cmd, '-sla DEADLINEsla')
    out_cmd = paste(out_cmd," \"",  cmd, "\"", sep = "")
    names(out_cmd)= names(cmd)
    return(out_cmd)
  }


##################
# Makes qsub command that wraps shell command "cmd" to send to queue "queue"
# redirebmccting output / error etc streams to path prefixed by "jname",
# optional_args: maximum memory requirements "mem", "jlabel" job label
##################
qsub_cmd = function(script.fn, queue, jname = NULL, jlabel = NULL, jgroup = NULL, mem=NULL, group = "cgafolk", cwd = NULL, mc.cores = NULL, deadline = F, now = FALSE)
    {
        if (is.null(jname) & is.null(names(script.fn)))
            jname = 'job'
        
        if (length(jname) != length(script.fn))
            jname = rep(jname, length(script.fn))
        
        if (!is.null(jname))
            names(script.fn) = dedup(jname)    
        
        qjname = paste( "\"", names(script.fn), "\"", sep="" )
        qjout = paste( "", names(script.fn), ".bsub.out", " " , sep="" )
        qjerr = paste( "", names(script.fn), ".bsub.err", "", sep="" )
        qjrout = paste( "", names(script.fn), ".R.out", "", sep="" )
                                        #        out_cmd = paste( "qsub -V -o ", qjout, " -e ",  qjerr)
        out_cmd = paste("qsub -V -j y -o ", qjout);
        out_cmd = paste(out_cmd, ifelse(is.na(queue), '', paste("-q ", queue)))
        ##        if (!is.null(mem)) out_cmd = paste(out_cmd, " -l mem=", mem, "gb", sep = "");
        if (!is.null(mem)) out_cmd = paste(out_cmd, " -l h_vmem=", mem, "g", sep = "");
        if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
        if (!is.null(cwd)) out_cmd = paste(out_cmd, " -wd ", cwd )
        if (!is.null(qjname)) out_cmd = paste(out_cmd, " -N ", jlabel)
        out_cmd = paste(out_cmd, '-now', ifelse(now, 'y', 'n'))
        if (!is.null(mc.cores)) out_cmd = paste(out_cmd, ifelse(!is.na(mc.cores), ifelse(mc.cores>1,  paste(" -pe smp",  mc.cores), ''), ''))
        out_cmd = paste(out_cmd, script.fn)
        names(out_cmd)= names(script.fn)
        return(out_cmd)
    }


#' @name Job-class
#' @rdname Job-class
#' @exportMethod report
#' @export
setGeneric('report', function(.Object, ...)  standardGeneric('report'))


#' @name report
#' @title Retrieves a data.table reporting detailed runtime stats associated with each job in the Job object.
#' @description
#'
#' Scrapes the underlying output directories associated with this Job object to determine whether jobs completed, stats on how much
#' memory and time was used by the job, when the jobs launched and when they finished, etc. 
#'
#' This is used by the Job update method, but can also be useful for debugging and diagnostics. 
#' 
#' @export
#' @author Marcin Imielinski
setMethod('report', 'Job', function(.Object, mc.cores = 1, force = FALSE)
    {
        out = cbind(data.table(runinfo(.Object)[, key(.Object), with = FALSE]), .parse.info(.Object@runinfo$stderr, mc.cores = mc.cores, force = force))
        suppressWarnings(out[ , key := NULL])
        setkeyv(out, key(.Object))
        return(out[1:nrow(out), ])
    })

.parse.info = function(jname, detailed = F, force = FALSE, mc.cores = 1)
{      

    
    dir = file.dir(jname)
    jname = file.name(jname)

  
  input.jname = jname
  jname = gsub('\\.bsub\\.out$', '', gsub('\\.bsub\\.err$', '', jname))
    names(input.jname) = jname
        
  if (length(jname)==0)    
    outs = data.frame(jname = NA,
      out.file = NA,
      err.file = NA,
      exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA,
        success = NA,
        stringsAsFactors = F)
  else
      {          
        outs = data.frame(jname = gsub('\\.R$', '', jname),
            out.file = paste(dir,'/', jname, '.bsub.out', sep = ''),
            err.file = paste(dir, '/', jname, '.bsub.err', sep = ''),
            exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA,
            success = NA,
            job_type = NA, 
            stringsAsFactors = F);

        fn = paste(dir, jname, '.bsub.out', sep = '')
        fn.err = paste(dir, jname, '.bsub.err', sep = '')
        fn.report = paste(dir, jname, '.bsub.report', sep = '')
        fn.report.sge = paste(dir, jname, '.bsub.report', sep = '')
        
        mtime = data.table(out = file.info(fn)$mtime, err = file.info(fn.err)$mtime, report = file.info(fn.report)$mtime, report.sge = file.info(fn.report.sge)$mtime)
        mtime[, report := pmax(report, report.sge, na.rm = TRUE)]

        ## we can use the report if the report exists and is younger than both the err and out
        ## or if (somehow) the err and out don't exist but the report does
        fn.rep.ex = mtime[ ,ifelse(!is.na(report), ifelse(!is.na(err) | is.na(out), pmin(report>err, report>out, na.rm = TRUE), FALSE), FALSE)] & !force
        
        if (any(fn.rep.ex))
            outs[fn.rep.ex, ] = do.call(rbind, lapply(fn.report[fn.rep.ex], read.delim, strings = FALSE))[, names(outs)]
        
        ## fn.ex these are the ones we need to parse again
        fn.ex = (file.exists(fn) | file.exists(fn.err)) & !fn.rep.ex; 

        if (!any(fn.ex))
            return(outs)
        
        tmp = matrix(unlist(mclapply(which(fn.ex),
            function(i)
                {
                    p = pipe(paste('tail -n 100', fn[i]))
                    y = readLines(p);
                    close(p)
                    p = pipe(paste('head -n 100', fn[i]))
                    sge = grep('FLOW', readLines(p), value = TRUE)
                    close(p)                    
                    if (any(grepl('^Sender.*LSF System', y))) ## LSF job
                        {
                            y = split(y, cumsum(grepl('^Sender', y)))
                            y = y[[length(y)]]  ## picks "last" dump from lsf to this out file
                            return(c('lsf',
                                     c(grep('^Exited with', y, value = T), grep('^Successfully completed', y, value = T), '')[1],
                                     c(grep('^TERM', y, value = T), '')[1],
                                     c(gsub('Started at ', '', grep('^Started at', y, value = T)), '')[1],
                                     c(gsub('Results reported ((at)|(on)) ', '', grep('^Results reported ((at)|(on))', y, value = T), ''))[1],
                                     c(gsub('[ ]+CPU time[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+CPU time', y, value = T)), '')[1],
                                     c(gsub('[ ]+Max Memory[ ]+\\:[ ]+(.*)', '\\1', grep('^[ ]+Max Memory', y, value = T)), '')[1],
                                     c(gsub('[ ]+Max Swap[ ]+\\:[ ]+(.*)', '\\1', grep('^[ ]+Max Swap', y, value = T)), '')[1],
                                        #                           c(gsub('[ ]+Max Memory[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+Max Memory', y, value = T)), '')[1],
                                        #                           c(gsub('[ ]+Max Swap[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+Max Swap', y, value = T)), '')[1],
                                     c(gsub('[ ]+Max Processes[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Processes', y, value = T)), '')[1],
                                     c(gsub('[ ]+Max Threads[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Threads', y, value = T)), '')[1]
                                     ))
                        }
                    else if (length(sge)>0)
                        {
                            fn.report.sge = paste(fn.report[i], '.sge', sep = '')
                            jobnum = gsub('^FLOW.SGE.JOBID=(.*)', '\\1',  sge[1])
                            p = pipe(paste('qacct -j', jobnum[length(jobnum)]))
                            tmp = readLines(p)
                            close(p)
                            vals = structure(str_trim(gsub('^\\S+\\s+(.*)', '\\1', tmp, perl = TRUE)), names = gsub('(^\\S+) .*', '\\1', tmp, perl = TRUE))
                            if (length(tmp)>0)                                
                                {
                                    write.table(as.data.frame(as.list(vals)), fn.report.sge, sep = '\t', quote = FALSE, row.names = FALSE)
                                }
                            else if (file.exists(fn.report.sge[i])) ## read from file if exists
                                {
                                    vals = unlist(read.delim(fn.report.sge, stringsAsFactors = FALSE))
                                }

                            cpuu= gsub('[^a-zA-Z]', '', vals['cpu'])
                            memu= gsub('[^a-zA-Z]', '', vals['mem'])
                            vmemu= gsub('[^a-zA-Z]', '', vals['maxvmem'])

                            return(c('sge',
                                     ifelse(vals['exit_status']=='0', 'Successfully completed.', vals['exit_status']),
                                     vals['failed'],
                                     vals['start_time'],
                                     vals['end_time'],
                                     as.numeric(gsub('[a-zA-Z]', '', vals['cpu']))*ifelse(grepl('[hH]', vmemu), 3600, ifelse(grepl('m', vmemu), 60, 1)),
                                     as.numeric(gsub('[a-zA-Z]', '',  vals['mem']))/ifelse(grepl('MB', vmemu), 1000, ifelse(grepl('KB', vmemu), 1e6, 1)),
                                     as.numeric(gsub('[a-zA-Z]', '', vals['maxvmem']))/ifelse(grepl('MB', vmemu), 1000, ifelse(grepl('KB', vmemu), 1e6, 1)),
                                     vals['slots'],
                                     vals['slots']
                                     ))
                        }
                    else ## interpret job as locally run with a /usr/bin/time -v output
                    {
                            y = tryCatch(readLines(fn.err[i]), error = function(e) NULL)
                            if (is.null(y))
                                y = readLines(fn[i])
                            ix = grep('Command being timed', y)
                            if (length(ix)==0) ## fail
                                return(rep(as.character(NA), 10))
                            ix = ix[length(ix)] ### only get the last instance                            
                            y = grep('\t.*', y[ix:length(y)][-1], value = TRUE)
                            tmp = strsplit(y, '[\\:\t]')
                            keyval = structure(str_trim(sapply(tmp, function(x) if (length(x)>2) x[[3]] else NA)),
                                names = sapply(tmp, function(x) if (length(x)>1) x[[2]] else NA))
                            etime = file.info(fn.err[i])$mtime
                            stime = etime - as.numeric(keyval['User time (seconds)'])
                            exit.status = ifelse(keyval['Exit status']==0, 'Successfully completed.', keyval['Exit status'])
                            return(c('local', exit.status, NA, as.character(stime), as.character(etime),
                                     keyval['User time (seconds)'], keyval['Maximum resident set size (kbytes)'], NA,
                                     as.numeric(gsub('\\%', '', keyval['Percent of CPU this job got']))/100, NA))
                        }
                }, mc.cores = mc.cores)), ncol = 10, byrow = T)


        
        colnames(tmp) = c('job.type', 'exit.flag', 'term.flag', 'started', 'reported', 'cpu.time', 'max.memory', 'max.swap', 'max.cpu', 'max.thr')
        
        .parse.mem = function(mem)
            {
                ix = !is.na(mem)
                out.mem = rep(NA, length(mem))
                if (any(ix))
                    {
                        mem = mem[ix]
                        tmp = strsplit(mem, '[ ]+')
                        tmp.mem = suppressWarnings(as.numeric(sapply(tmp, function(x) x[1])))
                        tmp.mem.units = sapply(tmp, function(x) x[2])
                        out.mem[ix] = ifelse(is.na(tmp.mem.units), tmp.mem/1e6/4, ## assume time output, which is in kbytes * 4
                                   ifelse(tmp.mem.units == 'MB', tmp.mem/1e3,
                                          ifelse(tmp.mem.units == 'KB', tmp.mem/1e6, 
                                                 tmp.mem)))

                    }
                return(out.mem)                
            }

        ## normalize to GB
        TIME.FORMAT1 = '%a %b %d %H:%M:%S %Y';
        TIME.FORMAT2 = '%Y-%m-%d %H:%M:%S';
        
        outs$job_type[fn.ex] = tmp[, 'job.type']
        outs$exit_flag[fn.ex] = tmp[, 'exit.flag']
        outs$term_flag[fn.ex] = tmp[, 'term.flag']
        outs$started[fn.ex] = ifelse( outs$job_type[fn.ex]%in% c('lsf','sge'),
                         as.character(as.POSIXct(strptime(tmp[, 'started'], TIME.FORMAT1))),
                        as.character(as.POSIXct(strptime(tmp[, 'started'], TIME.FORMAT2))))
         outs$reported[fn.ex] = ifelse(outs$job_type[fn.ex]%in% c('lsf','sge'),
                        as.character(as.POSIXct(strptime(tmp[, 'reported'], TIME.FORMAT1))),
                        as.character(as.POSIXct(strptime(tmp[, 'reported'], TIME.FORMAT2))))
        outs$hours_elapsed = as.numeric(as.POSIXct(outs$reported)-as.POSIXct(outs$started), units = 'hours')
        outs$cpu_time[fn.ex] = suppressWarnings(as.numeric(tmp[, 'cpu.time']))
        outs$max_mem[fn.ex] = ifelse(outs$job_type[fn.ex] == 'sge', as.numeric(tmp[, 'max.memory']), .parse.mem(tmp[, 'max.memory']))

        if (detailed)
            {
                outs$max_swap[fn.ex] = tmp[, 'max.swap']
                outs$max_processes[fn.ex] = tmp[, 'max.processes']      
                outs$max_threads[fn.ex] = tmp[, 'max.threads']
            }
        outs$success = ifelse(!is.na(outs$exit_flag), grepl('Success', outs$exit_flag), NA)
        rownames(outs) = dedup(outs$jname)

        ## cache row slices of this report table in the output directories
        ## for easier downstream access
        for (i in which(fn.ex))
            write.table(outs[i, ], fn.report[i], sep = '\t', quote = F, row.names = FALSE)
    }

    
  outs = as.data.table(outs)
  
  if (!is.null(input.jname))
      outs = outs[, key := input.jname[jname]]
  else
      outs = outs[, key := jname]

    setkey(outs, 'key')

  return(outs)  
}


#' @name xml2task 
#' @title Makes a best attempt to convert a firehose xml task configuration into a Task object or .task file
#' @description
#'
#' Takes a path to an xml file of a FH task configuration, module path, and outputs a Task object or writes
#' to a .task fil.e
#'
#' @exportMethod report
#' @export
#' @author Marcin Imielinski
xml2task = function(path, module = NULL, out.file = NULL)
    {
        require(XML)
        
        tasks = xmlToList(xmlParse(path))
        tasks = tasks[which(names(tasks)=='pipeline-configuration')]

        if (length(tasks)==0)
            stop('No pipeline configurations found in this xml file .. check file')
        
        out = lapply(tasks, function(task.config)
            {
                if (!is.null(module))           
                    if (!is(module, 'Module'))
                        module = Module(module)
                
                out = NULL
                if (!is.null(task.config$outputs))
                    if (length(task.config$outputs)>0)
                        {
                            outputs = as.data.table(do.call('rrbind', lapply(task.config$outputs, function(x) as.data.frame(rbind(unlist(x))))))        

                            setnames(outputs, gsub('\\-', '_', names(outputs)))

                            outs = outputs[, {           
                                list(list(FlowOutput(name = target_annotation_type_name, pattern= paste(expression, '.*', extension, "$", sep = ''))))
                            }, keyby = target_annotation_type_name]
                        }
                
                
                arg = NULL

                if (!is.null(task.config$parameter))
                    if (length(task.config$parameter)>0)
                        {
                            params = as.data.table(do.call('rrbind', lapply(task.config$"parameter", function(x) as.data.frame(rbind(unlist(x))))))
                            setnames(params, gsub('\\-', '_', names(params)))

                            if (is.null(params$"default_value"))
                                params$"default_value" = NA
                            
                            arg = params[, {
                                if (mode == 'FlowLiteral')
                                    list(list(FlowLiteral(name = name, arg = expression, path = file.exists(expression))))
                                else
                                    list(list(FlowAnnotation(name = name, arg = expression, path = file.exists(expression), default = default_value)))
                            }, keyby = name]
                        }

                if (is(module, 'Module'))
                    if (!is.null(arg))
                        if (!is.null(outs))
                            return(do.call(Task, c(structure(arg[[2]], names = arg[[1]]), list(outputs = structure(outs[[2]], names = outs[[1]]), module = module))))
                        else
                            return(do.call(Task, c(structure(arg[[2]], names = arg[[1]]), list(module = module))))
                    else
                        {
                            if (!is.null(outs))
                                return(do.call(Task, list(outputs = structure(outs[[2]]), names = outs[[1]], module = module)))
                            else
                                stop('No inputs or outputs in this task config, malformed task.config file?')
                        }
                else
                    {
                        warning('No module provided as input so just dumping mock task .task file to stdout')
                        out = paste('#', task.config$"name", task.config$"task-id")
                        out = c(out, '/path/to/module/directory')
                        
                        if (length(arg)>0)
                            out = c(out,
                                paste('input\t',
                                      arg[[1]],'\t', 
                                      ifelse(sapply(arg[[2]], is, 'FlowLiteral'), '"', ''),
                                      sapply(arg[[2]], function(x) x@arg),
                                      ifelse(sapply(arg[[2]], is, 'FlowLiteral'), '"', ''),
                                      '\t',
                                      ifelse(sapply(arg[[2]], function(x) x@path), 'path', 'value'),
                                      '\t',
                                      ifelse(sapply(sapply(arg[[2]], default), is.null), '', sapply(arg[[2]], default)),
                                      sep = ''))
                        
                        if (length(outs)>0)
                            out = c(out,
                                paste('output\t',
                                      sapply(outs[[2]], function(x) x@name), '\t',
                                      sapply(outs[[2]], function(x) x@pattern)))                


                        out = c(out, '')
                        if (is.null(out.file))                    
                            writeLines(out)
                        else
                            writeLines(out, out.file)

                        return(NULL)
                    }
            })

        if (all(sapply(out, is.null)))
            {
                cat('')
            }
        else
            {
                if (length(out)==1)
                    out = out[[1]]                
                return(out)
            }
    }

####################
## UTILITY FUNCTIONS
####################
# grabs filenames from list of paths
file.name = function(paths)
  {
    return(gsub('(^|(.*\\/))?([^\\/]*)', '\\3', paths))
  }

# grabs file.dirs from liOAst of paths
file.dir = function(paths)
  {
    return(gsub('(^|(.*\\/))?([^\\/]*)$', '\\2', paths))
  }

# relabels duplicates in a character vector with .1, .2, .3
# (where "." can be replaced by any user specified suffix)
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = setdiff(unique(x[dup]), NA)
  udup.ix = lapply(udup, function(y) which(x==y))
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)  
}



#' Merge Job output annotations with data.table
#'
#' @name merge
#' @aliases merge,Job-class
#' @rdname Job-class
#' @docType methods
#' @title Merges the output annotations associated with this job with another keyed data.table of entities.
#' @description
#'
#' As jobs complete, one may want to update a "master" data.table with the outputs of Jobs.  This can 
#' be useful for manual running of larger workflows to which a given task contributes. 
#' 
#' #@exportMethod merge
#' @param x data.table or Job
#' @param y data.table or Job
#' @param force logical flag whether to force overwrite
#' @param prefix prefix to add to columns merged from the Job
#' @param suffix suffix to add to columns merged from the Job
#' @param sep  separator to add to columns merged from the Job
#' @author Marcin Imielinski
#' @export
setGeneric('merge', function(x, y, ...) standardGeneric('merge'))
setMethod('merge', signature(x="Job", y = 'data.table'), function(x, y, suffix = NULL, prefix = NULL, force = FALSE, sep = '_') {
        if (!is.data.table(y))
            stop('y must be keyed data.table')

        if (is.null(data.table::key(y)))
            stop('y must be keyed data.table')
        
        if (data.table::key(y) != key(x))
            stop('y must be keyed data.table with same key as Job object')
        
        if (length(ov <- setdiff(intersect(names(outputs(x)), names(y)), key(x)))>0)
            col = unique(c(key(x), setdiff(names(y), ov)))
        else
            col = names(y)

        if (any(duplicated(names(y))))
            {
                warning('entities data.table has duplicate columns, deduping, check table')
                y = y[, unique(names(y)), with = FALSE]
            }

        if (any(duplicated(y[[key(x)]])))
            stop('Input table has duplicate instances of table key')
        
        ids = intersect(ids(x), y[[key(x)]])
        oids = setdiff(y[[key(x)]], ids)
        out = merge(y[, col, with = FALSE], outputs(x), by = key(x), all.x = TRUE)

        if (length(oids)>0) ## correct weird merge behavior in R
        {
            rn = dedup(out[[key(x)]])
            out = as.data.frame(out)
            rownames(out) = rn
            out[oids, colnames(y)] = as.data.frame(y[oids, colnames(y), with = FALSE])
            out = as.data.table(out)
            setkeyv(out, key(x))
        }
        

        if (length(ov)>0)
            {
                old = y[list(ids(x)), ov, with = FALSE]
                new = outputs(x)[, ov, with = FALSE]
                for (this.ov in ov)
                    {                        
                        ix = !is.na(new[[this.ov]]) | !is.na(old[[this.ov]])
                        ix[ix] = new[[this.ov]][ix] != old[[this.ov]][ix]
                        ix <- ifelse(is.na(ix), FALSE, ix)
                        if (any(ix))
                            {
                                old.mtime = file.info(old[[this.ov]][ix])$mtime
                                new.mtime = file.info(new[[this.ov]][ix])$mtime
                                ix2 <- ifelse(is.na(old.mtime>new.mtime), FALSE, old.mtime>new.mtime)
                                if (!force & any(ix2))
                                {
                                    warning('Newer annotations in external data.table are being over-written by new ones, keeping old annotations, call with force = TRUE to override')
                                    new[[this.ov]][ix][ix2] = old[[this.ov]][ix][ix2]
                                }

                                if (!force &
                                    any(ix <- is.na(new[[this.ov]]) & !is.na(old[[this.ov]]), na.rm = TRUE))
                                {
                                    warning('Existing annotations in external data.table are being over-written by NA annotations, keeping old annotations, call with force = TRUE to override')
                                    new[[this.ov]][ix2] = old[[this.ov]][ix2]
                                }
                            }
                        setkeyv(out, key(x))
                        out[ids(x),][[this.ov]] = new[[this.ov]]
                    }
            }

        ix = match(setdiff(names(outputs(x)), key(x)), names(out))
        if (!is.null(prefix))
            setnames(out, ix, paste(prefix, names(out)[ix], sep = sep))
        
        if (!is.null(suffix))
            setnames(out, ix, paste(names(out)[ix], suffix, sep = sep))

        setkeyv(out, data.table::key(y))
        return(out)
    })

setMethod('merge', signature(x='data.table', y="Job"), function(x, y, ...) {
        merge(y, x, ...)
    })
    


#' @name more
#' @title more
#'
#' @description
#' "more" +/- grep vector of files
#'
#' @param x vector of iles
#' @param grep string to grep in files (=NULL)
#' @author Marcin Imielinski
#' @export
more = function(x, grep = NULL, pipe = FALSE)
{
    if (is.null(grep))
        x = paste('more', paste(x, collapse = ' '))
    else
        x = paste('grep -H', grep, paste(x, collapse = ' '), ' | more')

    if (pipe)
        {
            p = pipe(x)
            out = readLines(p)
            close(p)
            return(out)
        }
    else
        system(x)
}

#' @name tailf
#' @title tailf
#'
#' @description
#' "tail -f" +/- grep vector of files
#'
#' @param x vector of iles
#' @param grep string to grep in files (=NULL)
#' @author Marcin Imielinski
#' @export
tailf = function(x, n = NULL, grep = NULL)
{
    if (is.null(grep))
        if (is.null(n))
            x = paste('tail -f', paste(x, collapse = ' '))
        else
            x = paste('tail -n', n, paste(x, collapse = ' '))
    else
        x = paste('grep -H', grep, paste(x, collapse = ' '), ' | more')
    system(x)
}



#############################
#' @name rrbind
#' @title rrbind
#'
#' @description
#'
#' like rbind, but takes the intersecting columns of the dfs
#'
#' if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa)
#' 
#' @param ... 
#' @author Marcin Imielinski
############################
rrbind = function(..., union = T)
  {     
    dfs = list(...);  # gets list of data frames
    if (any(ix <- sapply(dfs, function(x) class(x)[1])!='data.frame'))
        dfs[ix] = lapply(dfs[ix], as.data.frame)

    dfs = dfs[!sapply(dfs, is.null)]    
    dfs = dfs[sapply(dfs, ncol)>0]

    ## defactorize (need to do to cat without introducing NA's in weird places)
    dfs = lapply(dfs, function(x) { for (y in names(x)) if (is.factor(x[,y])) x[, y] = as.character(x[, y]); return(x)})
    
    names.list = lapply(dfs, names);
    classes = unlist(lapply(dfs, function(x) sapply(names(x), function(y) class(x[, y]))))
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    unshared.u = unique(unlist(unshared))
    ix = which(sapply(dfs, nrow)>0)
    expanded.dfs = lapply(ix, function(x)
      {
        dfs[[x]][, unshared[[x]]] = as.character(NA);
        return(dfs[[x]][, cols, drop = F])
      })
    
    out = do.call('rbind', expanded.dfs);
    
    if (any(uix <<- which(classes[unshared.u] != 'character')))
      {
          ix = match(unshared.u, names(out))
          for (j in uix) ### HACK to prevent stupid class mismatches leading to NA BS
              out[, ix[j]] = as(out[, ix[j]], classes[unshared.u[j]])
      }
    
    if (!union)
      {
        shared = setdiff(cols, unique(unlist(unshared)))
        out = out[, shared];
      }    
    
   return(out)
}


## #' @name merge
## #' @title Merges the output annotations associated with this job with another keyed data.table of entities.
## #' @description
## #'
## #' As jobs complete, one may want to update a "master" data.table with the outputs of Jobs.  This can 
## #' be useful for manual running of larger workflows to which a given task contributes. 
## #' 
## #' #@exportMethod merge
## #' @param x data.table or Job
## #' @param y data.table or Job
## #' @param force logical flag whether to force overwrite
## #' @param prefix prefix to add to columns merged from the Job
## #' @param suffix suffix to add to columns merged from the Job
## #' @param sep  separator to add to columns merged from the Job
## #' @author Marcin Imielinski
## #' @export

