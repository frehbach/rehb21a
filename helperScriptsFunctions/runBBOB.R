###
###
###
require(reticulate)

runCOCO <- function(solver,
                    instances = 1,
                    functions = 1:24,
                    dimensions=2,
                    number_of_batches=1,
                    current_batch=1, # 1, ..., number_of_batches
                    verbose=T,
                    solver_name="myalgorithm",
                    solverParameterList = list()
                    ){
    suite_name="bbob"
    solver_module=""
    max_runs=1e9
    budget = 100
    ## imports
    use_python("/usr/bin/python2", required = T)
    py_config()
    py_run_string("from __future__ import absolute_import, division, print_function")
    py_run_string("try: range = xrange
except NameError: pass")
    py_run_string("import os, sys")
    py_run_string("import time")
    py_run_string("import numpy as np")  # "pip install numpy" installs numpy"
    py_run_string("import cocoex")
    py_run_string("from scipy import optimize") # for tests with fmin_cobyla
    py_run_string("from cocoex import Suite, Observer, log_level")
    py_run_string("del absolute_import, division, print_function")
    py_run_string("try: import cma
except: pass")# cma.fmin is a solver option, "pip install cma" installs cma
    py_run_string("try: from scipy.optimize import fmin_slsqp
except: pass")# "pip install scipy" installs scipy
    py_run_string("try: range = xrange
except NameError: pass")# let range always be an iterator
    py_run_string("from cocoex import default_observers")  # see cocoex.__init__.py
    py_run_string("from cocoex.utilities import ObserverOptions, ShortInfo, ascetime, print_flush")
    py_run_string("from cocoex.solvers import random_search")


    ## default observer options
    py_run_string("def default_observer_options(budget_=None, suite_name_=None, current_batch_=None, solver_name_=None, solver_module_=None):
        global budget, suite_name, number_of_batches, current_batch
        if budget_ is None:
            budget_ = budget
        if suite_name_ is None:
            suite_name_ = suite_name
        if current_batch_ is None and number_of_batches > 1:
            current_batch_ = current_batch
        if solver_name_ is None:
            solver_name_ = solver_name
        if solver_module_ is None:
            solver_module_ = solver_module
        opts = {}
        try:
            opts.update({'result_folder': '\"%s_on_%s%s_budget%04dxD\"'
                        % (solver_name_,
                           suite_name_,
                           \"\" if current_batch_ is None
                              else \"_batch%03dof%d\" % (current_batch_, number_of_batches),
                           budget_)})
        except: pass
        try:
            opts.update({'algorithm_name': solver_name_ + solver_module_})
        except: pass
        return opts")
    ####


    ###### Parameters:
    py_run_string(paste("suite_name = '",suite_name,"'",sep=""))
    py_run_string(paste("budget = ",budget,sep=""))
    py_run_string(paste("max_runs = ",max_runs,sep=""))
    py_run_string(paste("number_of_batches = ",number_of_batches,sep=""))
    py_run_string(paste("current_batch = ",current_batch,sep=""))
    py_run_string(paste("solver_name = '",solver_name,"'",sep=""))
    py_run_string(paste("solver_module = '",solver_module,"'",sep=""))
    #py_run_string("suite_options = ''") #TODO: parameter
    py_run_string(paste("suite_options='dimensions: ",paste(dimensions,collapse = ","), " ",
                        "function_indices: ",paste(functions,collapse = ",")," ",
                        "instance_indices: ",paste(instances,collapse = ","),"'",sep=""))
    py_run_string("suite_instance = ''")
    py_run_string("observer_options = ObserverOptions({  # is (inherited from) a dictionary
        'algorithm_info': '\"A SIMPLE RANDOM SEARCH ALGORITHM\"'
    })")

    ## other stuff
    py_run_string("suite = Suite(suite_name, suite_instance, suite_options)")
    py_run_string("observer_name = default_observers()[suite_name]")
    py_run_string("observer_options.update_gracefully(default_observer_options())")
    py_run_string("observer = Observer(observer_name, observer_options.as_string)")
    ##batch loop:
    py_run_string("addressed_problems = []")
    py_run_string("short_info = ShortInfo()")
    py$budget
    ######
    py_run_string("problems = enumerate(suite)")
    py$suite$dimensions
    probs <- py$suite$problem_names
    np <- length(probs)
    if(verbose) print("number of problems:")
    if(verbose) print(np)
    batches <- sort(rep(1:number_of_batches,ceiling(np/number_of_batches)))
    for(i in 0:(np-1)){ #i is problem_index
        if(batches[i+1]!=current_batch){
            next
        }
        py_run_string(paste("problem = suite.get_problem(",i,")"))
        if(verbose) print(py$problem)
        py_run_string("observer.observe(problem)")
        #run optimizer with restarts??
        myfun <- py_to_r(py$problem)
        lower <- py$problem$lower_bounds
        upper <- py$problem$upper_bounds
        ## START!
        res <- solver(myfun,lower,upper,solverParameterList)
        py_run_string("problem.free()")
    }
}

startBatches <- function(nPerBatch,
                         solver,
                         instances = 1, dimensions=2,
                         verbose=T, solver_name="myalgorithm", solverParameterList){
    totalProblems <- capture.output(runCOCO(
        solver,current_batch = 2,number_of_batches = 1,
        dimensions=dimensions, instances = instances, solver_name = solver_name,solverParameterList))
    totalProblems <- as.numeric(stringr::str_split(totalProblems,stringr::fixed("[1] "))[[2]][2])
    nBatches <- ceiling(totalProblems/nPerBatch)
    for(i in 1:nBatches){
        runCOCO(
            solver,current_batch = i,number_of_batches = nBatches,solverParameterList = solverParameterList,
            dimensions=dimensions, instances = instances, solver_name = solver_name)
    }
}

