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

#' @name bsub_cmd
#' @title bsub_cmd
#' @description
#'
#' Makes bsub command that wraps shell command "cmd" to send to queue "queue"
#' redirecting output / error etc streams to path prefixed by "jname",
#' optional_args: maximum memory requirements "mem", "jlabel" job label
#'
#' @param cmd string shell command to be submitted via bsub
#' @param queue string name of specified queue, '-q "queue_name" ' (default = NULL)
#' @param jname  string path prefix by 'jname' (default = NULL)
#' @param jlabel string job name for '-J "job_name" ' (default = NULL)
#' @param jgroup string job_group_name for '-g "job_group_name" ' (default = NULL)
#' @param mem  integer amount of virtual memory/RAM to use via resource requirement arg '-R "res_req"'' (default = NULL)
#' @param group  string project_name for '-P' (default = NULL)
#' @param cwd string pathname to current working directory; -cwd "current_working_directory" (default = NULL)
#' @param mc.cores integer number of cores to use (default = 1)
#' @param deadline boolean specifies if deadline initiation time used (default = FALSE)
#' @author Marcin Imielinski
#' @export
bsub_cmd = function(cmd, queue = NULL, jname = NULL, jlabel = NULL, jgroup = NULL, mem = NULL, group = NULL, cwd = NULL, mc.cores = NULL, deadline = FALSE)
{
    if (is.null(jname) & is.null(names(cmd))){
        jname = 'job'
    }

    if (length(jname) != length(cmd)){
        jname = rep(jname, length(cmd))
    }

    if (!is.null(jname)){
        names(cmd) = dedup(jname)
    }

    qjname = paste( "\"", names(cmd), "\"", sep="" )
    qjout = paste( "\"", names(cmd), ".bsub.out", "\" ", sep="" )
    qjerr = paste( "\"", names(cmd), ".bsub.err", "\" ", sep="" )
    qjrout = paste( "\"", names(cmd), ".R.out", "\" ", sep="" )
    out_cmd = paste( "bsub -o ", qjout, " -e ",  qjerr);
    out_cmd = paste(out_cmd, ifelse(is.na(queue), '', paste("-q ", queue)))
#    if (!is.null(queue)) out_cmd = ifelse(is.na(queue), '', paste("-q ", queue))
    if (!is.null(group)) out_cmd = paste(out_cmd, " -P ", group)
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -R \"rusage[mem=", mem, "]\" ", sep = "")
    if (!is.null(jlabel)) out_cmd = paste(out_cmd, " -J ", jlabel )
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) out_cmd = paste(out_cmd, " -cwd ", cwd )
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, sprintf(" -n %d,%d -R 'span[hosts=1]'", mc.cores, mc.cores))
    if (deadline){ out_cmd = paste(out_cmd, '-sla DEADLINEsla') }
    out_cmd = paste(out_cmd," \"",  cmd, "\"", sep = "")
    names(out_cmd)= names(cmd)
    return(out_cmd)
}




#' @name qsub_cmd
#' @title qsub_cmd
#' @description
#'
#' Makes qsub command that wraps shell command "script.fn" to send to queue "queue"
#' redirecting output / error etc streams to path prefixed by "jname",
#' optional_args: maximum memory requirements "mem", "jlabel" job label
#'
#' @param script.fn string shell command to be submitted via bsub
#' @param queue string queue destination, specifies destination of the job '-q "destination" ' (default = NULL)
#' @param jname  string path prefix by 'jname' (default = NULL)
#' @param jlabel string job name for '-J "job_name" ' (default = NULL)
#' @param jgroup string job_group_name for '-g "job_group_name" ' (default = NULL)
#' @param mem  integer amount of virtual memory/RAM to use via resource requirement arg '-R "res_req"'' (default = NULL)
#' @param group  string project_name for '-P' (default = NULL)
#' @param cwd string pathname to current working directory; -cwd "current_working_directory" (default = NULL)
#' @param mc.cores integer number of cores to use (default = 1)
#' @param deadline boolean specifies if deadline initiation time used (default = FALSE)
#' @author Marcin Imielinski
#' @export
qsub_cmd = function(script.fn, queue = NULL, jname = NULL, jlabel = NULL, jgroup = NULL, mem = NULL, group = NULL, cwd = NULL, mc.cores = NULL, deadline = F, now = FALSE, touch_job_out = TRUE, qprior = 0)
{
    if (is.null(jname) & is.null(names(script.fn))){
        jname = 'job'
    }

    if (length(jname) != length(script.fn)){
        jname = rep(jname, length(script.fn))
    }

    if (!is.null(jname)){
        names(script.fn) = dedup(jname)
    }

    qjname = paste( "\"", names(script.fn), "\"", sep="" )
    qjout = paste( "", names(script.fn), ".bsub.out", " " , sep="" )
    qjerr = paste( "", names(script.fn), ".bsub.err", "", sep="" )
    qjrout = paste( "", names(script.fn), ".R.out", "", sep="" )
    out_cmd = paste0("qsub -V -j y -p ", qprior, " -o ", qjout)
    ## out_cmd = sprintf("qsub -V -j y -p %s -o %s ", qprior, qjout)
    ## out_cmd = paste("qsub -V -j y -o ", qjout);
    out_cmd = paste(out_cmd, ifelse(is.na(queue), '', paste("-q ", queue)))
    #if (!is.null(queue)) out_cmd = ifelse(is.na(queue), '', paste("-q ", queue))
    if (!is.null(group)) out_cmd = paste(out_cmd, " -P ", group)
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -l h_vmem=", mem, "g", sep = "")
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) {
        current_umask = Sys.umask(mode = NA)
        out_cmd = paste(out_cmd, " -wd ", cwd )
        Sys.umask(mode = "0002")
        base::file.create(trimws(paste0(cwd, "/", qjout)))
        Sys.umask(mode = current_umask)
        ## cmds = sprintf("umask 002; touch %s", paste0(cwd, "/", qjout))
        ## lapply(cmds, function(this_cmd) {system(this_cmd); return(NULL)})
    }
    if (!is.null(qjname)) out_cmd = paste(out_cmd, " -N ", jlabel)
    out_cmd = paste(out_cmd, '-now', ifelse(now, 'y', 'n'))
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, ifelse(!is.na(mc.cores), ifelse(mc.cores > 0,  paste(" -pe smp",  mc.cores), ''), ''))
    out_cmd = paste(out_cmd, script.fn)
    names(out_cmd)= names(script.fn)
    return(out_cmd)
}


#' @name ssub_cmd
#' @title ssub_cmd
#' @description
#'
#' Makes ssub command that wraps shell command "script.fn" to send to queue "queue"
#' redirecting output / error etc streams to path prefixed by "jname",
#' optional_args: maximum memory requirements "mem", "jlabel" job label
#'
#' @param script.fn string shell command to be submitted via bsub
#' @param queue string queue destination, specifies destination of the job '-q "destination" ' (default = NULL)
#' @param jname  string path prefix by 'jname' (default = NULL)
#' @param jlabel string job name for '-J "job_name" ' (default = NULL)
#' @param jgroup string job_group_name for '-g "job_group_name" ' (default = NULL)
#' @param mem  integer amount of virtual memory/RAM to use via resource requirement arg '-R "res_req"'' (default = NULL)
#' @param group  string project_name for '-P' (default = NULL)
#' @param cwd string pathname to current working directory; -cwd "current_working_directory" (default = NULL)
#' @param mc.cores integer number of cores to use (default = 1)
#' @param deadline boolean specifies if deadline initiation time used (default = FALSE)
#' @param now
#' @param time
#' @param qprior used for setting the "nice" value for SLURM
#' @param gres Sets gres value
#' @author Zoran Gajic
#' @export
ssub_cmd = function(script.fn, queue, qos = NULL, jname = NULL, jlabel = NULL, jgroup = NULL, mem=NULL, group = NULL, cwd = NULL, mc.cores = NULL, deadline = F, now = FALSE, time = "00", qprior = NULL, gres = NULL)
{
    if (is.null(jname) & is.null(names(script.fn)))
        jname = 'job'
    
    if (length(jname) != length(script.fn))
        jname = rep(jname, length(script.fn))

    if (is.null(qos))
        qos = as.character(NA)

    if (is.null(gres))
        gres = as.character(NA)
    
    if (!is.null(jname))
        names(script.fn) = dedup(jname)    
    qjname = paste( "\"", names(script.fn), "\"", sep="" )
    qjout = paste( "", names(script.fn), ".bsub.out", " " , sep="" )
    qjerr = paste( "", names(script.fn), ".bsub.err", "", sep="" )
    qjrout = paste( "", names(script.fn), ".R.out", "", sep="" )                    
    out_cmd = paste("sbatch --export=ALL --output=", qjout, " --error=", qjerr, sep = '');
    out_cmd = paste(out_cmd, ifelse(is.na(queue), '', paste0("-p ", queue)))
    out_cmd = paste(out_cmd, ifelse(is.na(qos), '', paste0("-q ", qos)))
    out_cmd = paste(out_cmd, ifelse(is.na(gres), "", paste0("--gres=", gres)))
    timestring = paste("--time=", time, ":00:00 ", sep = "")
    if (any(grepl(":", time))) timestring = paste("--time=", time, " ", sep = "")
    out_cmd = paste(out_cmd, timestring)
    if (!is.null(mem)) out_cmd = paste(out_cmd, " --mem=", mem, "G", sep = "");
    #if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    ## if (!is.null(cwd)) out_cmd = paste(out_cmd, " --workdir=", cwd , sep = '')
    if (!is.null(cwd)) {
        current_umask = Sys.umask(mode = NA)
        out_cmd = paste0(out_cmd, " --chdir=", cwd)
        Sys.umask(mode = "0002")
        base::file.create(trimws(paste0(cwd, "/", qjout)))
        Sys.umask(mode = current_umask)
        ## cmds = sprintf("umask 002; touch %s", paste0(cwd, "/", qjout))
        ## lapply(cmds, function(this_cmd) {system(this_cmd); return(NULL)})
    }
    if (!is.null(qjname)) out_cmd = paste(out_cmd, " --job-name=", jlabel, sep = '')
    ## out_cmd = paste(out_cmd, '-now', ifelse(now, 'y', 'n'))
    ## khadi Friday, Sep 25, 2020, Week 39, 08:59:16 AM
    ## wtf... slurm defaults to 2 cores per task if you don't set? changing >1 to >0
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, ifelse(!is.na(mc.cores), ifelse(mc.cores>0,  paste(" --cpus-per-task=",  mc.cores, sep = ""), ''), ''))
    if (!is.null(qprior)) out_cmd =  paste(out_cmd, ifelse(!is.na(qprior), ifelse(qprior>=0,  paste(" --nice=",  qprior, sep = ""), ''), '')) # set the "nice" value for SLURM
    out_cmd = paste(out_cmd, script.fn)
    names(out_cmd)= names(script.fn)
    return(out_cmd)
}



#' @name parse.info
#' @title parse.info
#' @description
#'
#' parses outputs of completed Job by subdirectory
#'
#' @param jname string 'jobname' pathanme to completed Job subdirectory
#' @param detailed boolean "verbose", outputs 'max.swap', 'max.processes', and 'max.threads' (default = FALSE)
#' @param force boolean force if *out/*err not found (default = FALSE)
#' @param mc.cores integer number of cores to use (default = 1)
#' @author Marcin Imielinski
#' @export
parse.info = function(jname, detailed = FALSE, force = FALSE, mc.cores = 1)
{
    dir = file.dir(jname)
    jname = file.name(jname)

    input.jname = jname
    jname = gsub('\\.bsub\\.out$', '', gsub('\\.bsub\\.err$', '', jname))
    names(input.jname) = jname

    if (length(jname)==0){
        outs = data.frame(jname = NA,
            out.file = NA,
            err.file = NA,
            exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA,
            success = NA,
            stringsAsFactors = F)
    }
    else{
        outs = data.frame(jname = gsub('\\.R$', '', jname),
            out.file = paste(dir,'/', jname, '.bsub.out', sep = ''),
            err.file = paste(dir, '/', jname, '.bsub.err', sep = ''),
            exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA,
            success = NA,
            job_type = NA,
            stringsAsFactors = F);
    }

    fn = paste(dir, jname, '.bsub.out', sep = '')
    fn.err = paste(dir, jname, '.bsub.err', sep = '')
    fn.report = paste(dir, jname, '.bsub.report', sep = '')
    fn.report.sge = paste(dir, jname, '.bsub.report', sep = '')

    mtime = data.table(out = file.info(fn)$mtime, err = file.info(fn.err)$mtime, report = file.info(fn.report)$mtime, report.sge = file.info(fn.report.sge)$mtime)
    mtime[, report := pmax(report, report.sge, na.rm = TRUE)]

    ## we can use the report if the report exists and is younger than both the err and out
    ## or if (somehow) the err and out don't exist but the report does
    fn.rep.ex = mtime[ ,ifelse(!is.na(report), ifelse(!is.na(err) | is.na(out), pmin(report>err, report>out, na.rm = TRUE), FALSE), FALSE)] & !force

    if (any(fn.rep.ex)){
        outs[fn.rep.ex, ] = do.call(rbind, lapply(fn.report[fn.rep.ex], read.delim, strings = FALSE))[, names(outs)]
    }

    ## fn.ex these are the ones we need to parse again
    fn.ex = (file.exists(fn) | file.exists(fn.err)) & !fn.rep.ex;

    if (!any(fn.ex)){
        return(outs)
    }

    tmp = matrix(unlist(mclapply(which(fn.ex),
        function(i){
            p = pipe(paste('tail -n 100', fn[i]))
            y = readLines(p);
            close(p)
            p = pipe(paste('head -n 100', fn[i]))
            sge = grep('FLOW', readLines(p), value = TRUE)
            close(p)
            if (any(grepl('^Sender.*LSF System', y))){   ## LSF job
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
                    c(gsub('[ ]+Max Processes[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Processes', y, value = T)), '')[1],
                    c(gsub('[ ]+Max Threads[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Threads', y, value = T)), '')[1]
                ))
            }
            else if (length(sge)>0){
                fn.report.sge = paste(fn.report[i], '.sge', sep = '')
                jobnum = gsub('^FLOW.SGE.JOBID=(.*)', '\\1',  sge[1])
                p = pipe(paste('qacct -j', jobnum[length(jobnum)]))
                tmp = readLines(p)
                close(p)
                vals = structure(str_trim(gsub('^\\S+\\s+(.*)', '\\1', tmp, perl = TRUE)), names = gsub('(^\\S+) .*', '\\1', tmp, perl = TRUE))
                if (length(tmp)>0){
                    write.table(as.data.frame(as.list(vals)), fn.report.sge, sep = '\t', quote = FALSE, row.names = FALSE)
                }
                ## read from file if exists
                else if (file.exists(fn.report.sge[i])){
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
            ## interpret job as locally run with a /usr/bin/time -v output
            else{
                y = tryCatch(readLines(fn.err[i]), error = function(e) NULL)
                if (is.null(y)){
                    y = readLines(fn[i])
                }
                ix = grep('Command being timed', y)
                if (length(ix)==0){
                    ## fail
                    return(rep(as.character(NA), 10))
                }
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

    .parse.mem = function(mem){
        ix = !is.na(mem)
        out.mem = rep(NA, length(mem))
        if (any(ix)){
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

    if (detailed){
        outs$max_swap[fn.ex] = tmp[, 'max.swap']
        outs$max_processes[fn.ex] = tmp[, 'max.processes']
        outs$max_threads[fn.ex] = tmp[, 'max.threads']
    }

    outs$success = ifelse(!is.na(outs$exit_flag), grepl('Success', outs$exit_flag), NA)
    rownames(outs) = dedup(outs$jname)

    ## cache row slices of this report table in the output directories
    ## for easier downstream access
    for (i in which(fn.ex)){
        write.table(outs[i, ], fn.report[i], sep = '\t', quote = F, row.names = FALSE)
    }

    outs = as.data.table(outs)

    if (!is.null(input.jname)){
      outs = outs[, key := input.jname[jname]]
    }
    else{
      outs = outs[, key := jname]
    }

    setkey(outs, 'key')

    return(outs)
}




#' @name xml2task
#' @title Makes a best attempt to convert a firehose xml task configuration into a Task object or .task file
#' @description
#'
#' Takes a path to an xml file of a FH task configuration, module path, and outputs a Task object or writes
#' to a .task file
#'
#' @param path string of pathname to firehose xml task
#' @param module string name of Module (default = NULL)
#' @param out.file string pathname to *.task file (default = NULL)
#' @author Marcin Imielinski
#' @export
xml2task = function(path, module = NULL, out.file = NULL)
{
    if (!requireNamespace("XML", quietly = TRUE)) {
        stop('In order to read xml files you must install the package "XML".')
    }

    tasks = XML::xmlToList(XML::xmlParse(path))
    tasks = tasks[which(names(tasks)=='pipeline-configuration')]

    if (length(tasks)==0){
        stop('No pipeline configurations found in this xml file .. check file')
    }

    out = lapply(tasks, function(task.config){
        if (!is.null(module)){
            if (!is(module, 'Module')){
                module = Module(module)
            }
        }

        out = NULL
        if (!is.null(task.config$outputs)){
            if (length(task.config$outputs)>0){
                outputs = as.data.table(do.call('rrbind', lapply(task.config$outputs, function(x) as.data.frame(rbind(unlist(x))))))

                setnames(outputs, gsub('\\-', '_', names(outputs)))

                outs = outputs[, {list(list(FlowOutput(name = target_annotation_type_name, pattern= paste(expression, '.*', extension, "$", sep = ''))))},
                    keyby = target_annotation_type_name]
            }
        }

        arg = NULL

        if (!is.null(task.config$parameter)){
            if (length(task.config$parameter)>0){
                params = as.data.table(do.call('rrbind', lapply(task.config$"parameter", function(x) as.data.frame(rbind(unlist(x))))))
                setnames(params, gsub('\\-', '_', names(params)))

                if (is.null(params$"default_value")){
                    params$"default_value" = NA
                }

                arg = params[, {if (mode == 'FlowLiteral'){
                                    list(list(FlowLiteral(name = name, arg = expression, path = file.exists(expression))))
                                }
                                else{
                                    list(list(FlowAnnotation(name = name, arg = expression, path = file.exists(expression), default = default_value)))
                                }
                            }, keyby = name]
            }
        }

        if (is(module, 'Module')){
            if (!is.null(arg)){
                if (!is.null(outs)){
                    return(do.call(Task, c(structure(arg[[2]], names = arg[[1]]), list(outputs = structure(outs[[2]], names = outs[[1]]), module = module))))
                }
                else{
                    return(do.call(Task, c(structure(arg[[2]], names = arg[[1]]), list(module = module))))
                }
            }
            else{
                if (!is.null(outs)){
                    return(do.call(Task, list(outputs = structure(outs[[2]]), names = outs[[1]], module = module)))
                }
                else{
                    stop('No inputs or outputs in this task config, malformed task.config file?')
                }
            }
        }
        else{
            warning('Warning: No module provided as input so just dumping mock task .task file to stdout')
            out = paste('#', task.config$"name", task.config$"task-id")
            out = c(out, '/path/to/module/directory')

            if (length(arg)>0){
                out = c(out, paste('input\t', arg[[1]],'\t',
                    ifelse(sapply(arg[[2]], is, 'FlowLiteral'), '"', ''),
                    sapply(arg[[2]], function(x) x@arg),
                    ifelse(sapply(arg[[2]], is, 'FlowLiteral'), '"', ''),
                    '\t',
                    ifelse(sapply(arg[[2]], function(x) x@path), 'path', 'value'),
                    '\t',
                    ifelse(sapply(sapply(arg[[2]], default), is.null), '', sapply(arg[[2]], default)),
                    sep = ''))
            }
            if (length(outs)>0){
                out = c(out,paste('output\t',
                    sapply(outs[[2]], function(x) x@name), '\t',
                    sapply(outs[[2]], function(x) x@pattern)))
            }

            out = c(out, '')
            if (is.null(out.file)){
                writeLines(out)
            }
            else{
                writeLines(out, out.file)
            }

            return(NULL)
        }
    })

    if (all(sapply(out, is.null))){
        cat('')
    }
    else{
        if (length(out)==1){
            out = out[[1]]
        }
        return(out)
    }

}




#' @name file.name
#' @title file.name
#' @description
#'
#' grabs filenames from list of paths
#'
#' @param paths vector of pathnames
#' @author Marcin Imielinski
#' @export
file.name = function(paths)
{
    return(gsub('(^|(.*\\/))?([^\\/]*)', '\\3', paths))
}




#' @name file.dir
#' @title file.dir
#' @description
#'
#' grabs file.dirs from list of paths
#'
#' @param paths vector of pathnames
#' @author Marcin Imielinski
#' @export
file.dir = function(paths)
{
    return(gsub('(^|(.*\\/))?([^\\/]*)$', '\\2', paths))
}




#' @name dedup
#' @title dedup
#' @description
#'
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x vector of files
#' @param suffix string suffix (default = '.')
#' @author Marcin Imielinski
#' @export
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




#' @name more
#' @title more
#'
#' @description
#' "more" +/- grep vector of files
#'
#' @param x vector of files
#' @param grep string to grep in files (default = NULL)
#' @param pipe boolean grep piped input 'x' (default = NULL)
#' @author Marcin Imielinski
#' @export
more = function(x, grep = NULL, pipe = FALSE)
{
    if (is.null(grep)){
        x = paste('more', paste(x, collapse = ' '))
    }
    else{
        x = paste('grep -H', grep, paste(x, collapse = ' '), ' | more')
    }
    if (pipe){
        p = pipe(x)
        out = readLines(p)
        close(p)
        return(out)
    }
    else{
        system(x)
    }
}




#' @name tailf
#' @title tailf
#'
#' @description
#' "tail -f" +/- grep vector of files
#'
#' @param x vector of files
#' @param n integer specific number of lines in file (default = NULL)
#' @param grep string to grep in files (default = NULL)
#' @author Marcin Imielinski
#' @export
tailf = function(x, n = NULL, grep = NULL)
{
    ## create command
    if (is.null(grep)){
        if (is.null(n)){
            x = paste('tail -f', paste(x, collapse = ' '))
        }
        else{
            x = paste('tail -n', n, paste(x, collapse = ' '))
        }
    }
    else{
      x = paste('grep -H', grep, paste(x, collapse = ' '), ' | more')
    }

    ## execute command
    system(x)
}


deparse1 = function (expr, collapse = " ", width.cutoff = 500L, ...) {
    base::paste(base::deparse(expr, width.cutoff, ...), collapse = collapse)
}


#' get S4 Slot
#' 
#' Simple function that returns NULL if slot not present
#' instead of erroring out
#' 
#' @export
getslot = function(object, name, default = NULL) {
    suppressWarnings({
        is_character = tryCatch(is.character(name), error = function(e) FALSE)
        is_symbol_name = tryCatch(is.symbol(name), error = function(e) FALSE)
        is_symbol_substitute_name = tryCatch(is.symbol(substitute(name)), error = function(e) FALSE)
        is_to_be_deparsed = is_symbol_name || is_symbol_substitute_name
    })
    if (is_character) {
        nm = name
    }
    if (is_to_be_deparsed) {
        if (is_symbol_substitute_name) nm = substitute(name)
        nm = deparse1(nm)
    }
    ## Note object can be NULL
    if (!methods::.hasSlot(object, nm)) return(default)
    return(methods::slot(object, nm))
}

#' Get chain of S4 methods
#' 
#' Function that goes through the chain of arguments and returns a default
#' instead of erroring out
#' 
#' @export
getslotchain = function(object, ..., default = NULL) {
    args = match.call(expand.dots = FALSE)$...
    out_object = object
    for (arg in args) {
        out_object = getslot(out_object, as.character(arg), default = default)
    }
    return(out_object)
}

`%@%` = getslot


#' @name reset.job
#' @title reset.job
#'
#' Reset a job with different params
#'
#' @return A Flow job object
#' @author Kevin Hadi
#' @export reset.job
reset.job = function(x, delete_cache = FALSE, ..., i = NULL, rootdir, jb.mem, jb.cores, jb.time, update_cores = 1, shell, force_shell, force_profile, task = NULL) {
    if (missing(jb.time))
        jb.time = base::get0("time", as.environment(x@runinfo), ifnotfound = "24:00:00")
    
    if (missing(jb.mem))
        jb.mem = base::get0("mem", as.environment(x@runinfo), ifnotfound = 16)
    
    if (missing(jb.cores))
        jb.cores = base::get0("cores", as.environment(x@runinfo), ifnotfound = 1L)

    if (missing(rootdir)) 
        rootdir = x@rootdir
    
    if (missing(force_shell)) {
        force_shell = Flow::getslotchain(x, task, module, force_shell, default = TRUE)
    }

    if (missing(force_profile)) {
        force_profile = Flow::getslotchain(x, task, module, force_profile, default = TRUE)
    }
    
    if (missing(shell)) {
        shell = Flow::getslotchain(x, task, module, shell, default = "bash")
    }
    
    if (!inherits(x, "Job")) stop ("x must be a Flow Job object")

    if (is.null(task)) {
        usetask = x@task
        ## empty_profile_slot = list(FlowAnnotation(name = "profile", arg = list(), path = FALSE, default = NA_character_))
        empty_profile_slot = NULL
        ## If profiles is not present - it's NULL
        profiles = getslotchain(usetask, profiles, default = empty_profile_slot)
        usetask@profiles = profiles
        usetask@module@shell = shell
        ## usetask@module@force_shell = getslotchain(usetask, module, shell, default = force_shell)
        usetask@module@force_shell = force_shell
        usetask@module@force_profile = force_profile
    } else if (is.character(task) || inherits(task, "Task"))
        usetask = task
    
    args = list(...)
    new.ent = copy(entities(x))
    
    if (!is.null(i)) {
        jb.mem = replace(x@runinfo$mem, i, jb.mem)
        jb.cores = replace(x@runinfo$cores, i, jb.cores)
    }

    tbl_task = viewtask(usetask)
    ## if (!all(names(args) %in% colnames(new.ent)))
    if (!all(names(args) %in% colnames(new.ent)) && !names(args) %in% tbl_task$V2)
        stop("adding additional column to entities... this function is just for resetting with new arguments")
    
    for (j in seq_along(args)) {
        data.table::set(new.ent, i = i, j = names(args)[j], value = args[[j]])
    }
    these.forms = formals(body(findMethods("initialize")$Job@.Data)[[2]][[3]])

    path_to_cache = Flow::getcache(x)

    if (delete_cache) {
        message("DELETING PATH TO .RDS:")
        message(path_to_cache)
        system2("rm", c("-f", path_to_cache))
    }

    jb = Job(
        task = usetask, 
        entities = new.ent, 
        rootdir = rootdir, 
        mem = jb.mem, 
        time = jb.time, 
        cores = jb.cores, 
        update_cores = update_cores,
        shell = shell,
        force_shell = force_shell,
        force_profile = force_profile
    )

    # if ("time" %in% names(these.forms)) {
    #     if ("update_cores" %in% names(these.forms))
    #         jb = Job(usetask, new.ent, rootdir = rootdir, mem = jb.mem, time = jb.time, cores = jb.cores, update_cores = update_cores)
    #     else
    #         jb = Job(usetask, new.ent, rootdir = rootdir, mem = jb.mem, time = jb.time, cores = jb.cores)
    # } else {
    #     if ("update_cores" %in% names(these.forms))
    #         jb = Job(usetask, new.ent, rootdir = rootdir, mem = jb.mem, cores = jb.cores, update_cores = update_cores)
    #     else
    #         jb = Job(usetask, new.ent, rootdir = rootdir, mem = jb.mem, cores = jb.cores)
    # }
    return(jb)
}


#' @name silent
#' @title run expression without any printed output
#'
#' execute expression without any output to console.
#' silent({var = function_that_has_explicit_print(...)})
#' 
#'
#' @author Kevin Hadi
#' @param ... an expression
#' @return NULL
#' @export
silent = function(this_expr, this_env = parent.frame()) {
    eval(expr = {capture.output(
            capture.output(... = this_expr,
                           file = "/dev/null",
                           type = c("output")),
            file = "/dev/null",
            type = "message")
    }, envir = this_env)
    invisible()
}

#' or Flow Task object
#'
#'
#' @author Kevin Hadi
#' @return A Flow job object
#' @export
viewtask = function(jb, arglst = c("name", "arg", "default")) {
    lst.emptychar2na = function(x) {
        x[!nzchar(x)] = NA_character_
        x
    }
    lst.zerochar2empty = function(x) {
        x[x == "character(0)"] = list("")
        x
    }
    ifun = function(x, arglst = arglst) {
        
        unlist(lst.emptychar2na(lst.zerochar2empty(lapply(arglst, function(y)
            tryCatch((slot(x, y)), error = function(e) NA_character_)))))
    }
    if (inherits(jb, "Job"))
        obj = jb@task
    else if (inherits(jb, "Task"))
        obj = jb
    else if (inherits(jb, "character"))
        obj = Task(jb)
    as.data.table(data.table::transpose(lapply(obj@args, ifun, arglst = arglst)))
}


#' getcache
#'
#' Get the path of a Flow job object's cache.
#'
#' @return A character
#' @export
getcache = function(object) {
      path = paste(object@rootdir, "/", task(object)@name,
                   ".rds", sep = "")
      return(path)
}
