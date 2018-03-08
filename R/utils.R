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
##    if (!is.null(queue)) out_cmd = ifelse(is.na(queue), '', paste("-q ", queue)) 
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
qsub_cmd = function(script.fn, queue = NULL, jname = NULL, jlabel = NULL, jgroup = NULL, mem = NULL, group = NULL, cwd = NULL, mc.cores = NULL, deadline = F, now = FALSE)
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
    out_cmd = paste("qsub -V -j y -o ", qjout);
    out_cmd = paste(out_cmd, ifelse(is.na(queue), '', paste("-q ", queue)))
    #if (!is.null(queue)) out_cmd = ifelse(is.na(queue), '', paste("-q ", queue))
    if (!is.null(group)) out_cmd = paste(out_cmd, " -P ", group)
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -l h_vmem=", mem, "g", sep = "")
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) out_cmd = paste(out_cmd, " -wd ", cwd )
    if (!is.null(qjname)) out_cmd = paste(out_cmd, " -N ", jlabel)
    out_cmd = paste(out_cmd, '-now', ifelse(now, 'y', 'n'))
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, ifelse(!is.na(mc.cores), ifelse(mc.cores > 1,  paste(" -pe smp",  mc.cores), ''), ''))
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
    require(XML)
        
    tasks = xmlToList(xmlParse(path))
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
rrbind = function(..., union = TRUE){     
    dfs = list(...);  # gets list of data frames
    if (any(ix <- sapply(dfs, function(x) class(x)[1])!='data.frame')){
        dfs[ix] = lapply(dfs[ix], as.data.frame)
    }

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
    expanded.dfs = lapply(ix, function(x){
        dfs[[x]][, unshared[[x]]] = as.character(NA);
        return(dfs[[x]][, cols, drop = F])
    })
    
    out = do.call('rbind', expanded.dfs);
    
    if (any(uix <<- which(classes[unshared.u] != 'character'))){
        ix = match(unshared.u, names(out))
        for (j in uix){
            ### HACK to prevent stupid class mismatches leading to NA BS
            out[, ix[j]] = as(out[, ix[j]], classes[unshared.u[j]])
        }
    }
    
    if (!union){
        shared = setdiff(cols, unique(unlist(unshared)))
        out = out[, shared];
    }    
    
   return(out)

}



