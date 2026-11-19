#! /usr/bin/env python3
"""
Usage: python3 imsil.py [options] [input-file ...]

The purpose of this script is to run imsil, in the following referred to as PROGRAM

options    : -v <version>       : Version of PROGRAM to be used (name of the subdirectory of 
                                  the PROGRAM code directory). Default is where the current
                                  directory is subdirectory of
             -p <parallel_jobs> : Maximum number of jobs (input files) to run in parallel
                                  parallel_jobs=0 means run all jobs in parallel
                                  (default is parallel_jobs=1)
             -t <threads>       : Number of OpenMP threads per job 
                                  threads=0 means use sequential version
                                  (default is threads=0)
             -r <hostname>      : host to run PROGRAM on
             -n <niceness>      : nice command for Unix. <niceness> means the niceness value 
                                  (generally between 1 and 20)
             -l <logfile>       : logfile, replacing standard output or nohup.out
             -s                 : silent execution (no questions asked)
             -i                 : interactive (foreground) execution
             -d                 : debug
             -D                 : use code with debug information (default if -d is set)
             
             If logfile exists, output is appended.
                
             Unless -s is set, the user is prompted for modifying these options.
 
input-file : the name of an input files.  The input file must have the extension '.inp' or 
             must be named 'INP'. '.inp' or 'INP' can be omitted on the command line. The 
             path of the input file must either be specified as an absolute path or relative 
             to the directory where the script is started.

If more than one input file is specified, PROGRAM is run with each input file, i.e., a "job"
is created for each file. Each job can use the given number of threads (-t option), resulting 
in a number of processes equal to the number of jobs times the number of threads per job.   
             
If the input file is given as an '.inp' file (as opposed to 'INP'), the input file and 
all files with the same name but different extension are copied into a temporary working 
directory with the name equal to the name of the files without extension (inflate). The 
names of the files in the new directory are identical to their original extensions 
converted to uppercase letters. After completion of PROGRAM, the reverse copy process 
(deflate) is executed on all all-capital files in the directory including the files generated 
by PROGRAM; if no not-all-capital files and no subdirectories have been created, the directory
is deleted. In this way, parallel runs of PROGRAM on input files in the same directory do 
not interfere. As an example
 
    input file name: test.inp
    start PROGRAM in interactive mode: python3 imsil.py test
    temporary working directory: test
    copy of input file: test/INP
    PROGRAM generates file OUT => test/OUT
    script copies file back to original directory as: test.out
 
In addition, a file PATH is written to the working directory containing the full paths
of directories where PROGRAM may look for additional input files:
   
&path
inppath='***/'         (path of original input file)
workpath='***/'        (working directory = path of input file copy INP)
tablepath='***/'       (PROGRAM root directory/tables)
/

It may be read by a FORTRAN namelist read statement.
""" 

import argparse, functools, glob, multiprocessing, runpy, os, shutil, socket, \
       subprocess, sys, time

PROGRAM = os.path.splitext( os.path.basename( __file__ ) )[0]

if os.name == 'posix':
    installdir = os.path.join( os.path.expanduser('~'), 'work', PROGRAM )  
                                                # parent directory of code and tables,
                                                # which contain the executable and the tables,
                                                # respectively
    #debugcmd = '/usr/bin/ddd' # full path of debug command
    debugcmd = '/usr/bin/nemiver' # full path of debug command
    #debugcmd = '/opt/intel/compilers_and_libraries_2016.3.210/linux/bin/intel64/gdb-ia' # full path of debug command
    #debugcmd = '/opt/intel/composer_xe_2013_sp1.3.174/bin/intel64/idb' # full path of debug command
else:                                   # unused for a long time
    installdir = 'e:\\Programme'
    debugcmd = r'C:\Program Files\Microsoft Visual Studio\Common\MSDEV98\BIN\dfdev.exe'

#############################################################################################
def program_string( version ):
    """
    Return program and version information as a string.
    Meant for use in output messages.
    """
    if os.path.isdir( os.path.join( installdir, version ) ):
        return PROGRAM.upper( ) + '(' + version + ')'
    else:
        return PROGRAM.upper( )

#############################################################################################

def parseCommandLine( args, *exts ):
    """
    Parse the command line.
       
    Arguments : args : list of command line arguments
                exts : list of allowed extensions (may also specify file in subdirectory
                       like '/INP')
       
    Results   : options : command line options, including inpFiles = list of input files 
                          (absolute or relative paths, with or without extension)
    """  
    
    # define the options
    
    # default host is the local host
    current_host = socket.gethostname()
    
    p = argparse.ArgumentParser( )
    p.add_argument( 'inpFiles', nargs='*' )
    p.add_argument( '--version', '-v' )
    p.add_argument( '--niceness', '-n', type=int, default=0 )
    p.add_argument( '--logfile', '-l', default=None )
    p.add_argument( '--silent', '-s', action='store_true' )
    p.add_argument( '--background', action='store_false' )
    p.add_argument( '--nobackground', '-i', action='store_true' )
    p.add_argument( '--debug', '-d', action='store_true' )
    p.add_argument( '--nodebug', action='store_true' )
    p.add_argument( '--debugcode', '-D', action='store_true' )
    p.add_argument( '--nodebugcode', action='store_true' )
    p.add_argument( '--parallel_jobs', '-p', type=int, default=1 )
    p.add_argument( '--threads', '-t', type=int, default=0 )
    p.add_argument( '--remote', '-r', dest='remotehost', default=current_host )
    p.add_argument( '--local', dest='localhost', default=current_host )
    p.add_argument( '--startdir', default=os.getcwd( ) )
    
    # parse the options
    
    options = p.parse_args( args[1:] )

    # default version (cannot be done as an argument to add_argument since it needs 
    # options.startdir)

    if not options.version:
        relcwd = os.path.relpath( options.startdir, installdir )
        options.version = relcwd.split(os.sep)[0]
    
    # override choices by '--no...' options
    
    if options.nobackground: options.background=False
    if options.nodebug: options.debug=False
    if options.nodebugcode: options.debugcode=False
    
    # prevent impossible options on Windows machines
    
    if os.name != 'posix':
        if options.background:
            options.background = False
            print( 'Running in background is only possible under Linux. ' +
                   'Running in foreground instead...' )
        if options.remote_host != current_host:
            options.remote_host = current_host
            print( 'Running on a remote host is only possible under Linux. ' +
                   'Running locally instead...' )
        if options.niceness:
            options.niceness = 0
            print( 'Niceness only possible under Linux. Ignoring niceness...' )

    # if no input files have been specified, ask for them (<CR> ends input)
 
    if not options.inpFiles and not options.silent:
        prompt = '\nPlease specify the names of the input files (' + exts[0] 
        for ext in exts[1:]:
            prompt = prompt + ' or ' + ext
        prompt = prompt + '). <CR> ends input.' 
        print( prompt )
        while True:
            inpFile = input( 'Input file: ' )
            if not inpFile: 
                break
            options.inpFiles.append( inpFile )
 
    if not options.inpFiles:
        sys.exit( args[0] + ': No input files given.' )
    
    # ask if other options should be modified

    if not options.silent:
        while 1:
            prompt = '\nCurrent options are:\n'
            if options.background:
                prompt = prompt + '\n   --background   ... background execution'
            else:
                prompt = prompt + '\n   -i             ... interactive (foreground) execution' 
            if options.debug:
                prompt = prompt + '\n   -d             ... run with debugger'
            else:
                prompt = prompt + '\n   --nodebug      ... run without debugger '
            if options.debugcode:
                prompt = prompt + '\n   -D             ... use code with debug information'
            else:
                prompt = prompt + '\n   --nodebugcode  ... use optimized code'
            prompt = prompt + '\n\nEnter "i", "d", or "D" to invert these options,'
            prompt = prompt + ' or <CR> to accept: '
            answer = input( prompt )
            if not answer: break
            if answer == 'i':
                options.background = not options.background
            elif answer == 'd':
                options.debug = not options.debug
                options.debugcode = not options.debugcode
            elif answer == 'D':
                options.debugcode = not options.debugcode
    
    return options

#############################################################################################

def setProgramPaths( options ):
    """
    Set various paths and check on existence.
       
    Arguments : options : command line options
  
    Results   : exeFile   : full path of PROGRAM executable
                tablesDir : full path of directory containing additional input files
    
    The code directory is PROGRAM/version/code or, if it does not exist, PROGRAM/code.
    The tables directory must be a sibling of the code directory
    """
    
    # code directory
    
    codedir = os.path.join( installdir, options.version, 'code' )
    
    if not os.path.isdir( codedir ):
        codedir_alt = os.path.join( installdir, 'code' )   # alternate code directory
        if os.path.isdir( codedir_alt ):
            codedir = codedir_alt
        else:
            sys.exit( program_string( options.version ) 
                      + ' possible code directories do not exist: ' 
                      + '\n   ' + codedir 
                      + '\n   ' + codedir_alt
                      + '\nSpecify valid version, adjust path in start script, ' 
                      + '\nor install code in different place' )
    
    # tables directory
    
    tablesDir = os.path.join( codedir, '..', 'tables' )
    if not os.path.isdir( tablesDir ):
        print( 'Warning: Tables directory does not exist:' + tablesDir )
    
    # path of executable file
    
    exeFile = 'imsil'
    if options.debugcode:
        exeFile = exeFile + '_db'
    if options.threads != 0:
        exeFile = exeFile + '_omp'
    exeFile = os.path.join( codedir, exeFile + '.exe' )
    
    if os.name != 'posix' and options.debug and options.debugcode:
        sys.exit( "exeFile = codedir/imsil_db.dsw?" )
    
    if not os.path.isfile( exeFile ):
        sys.exit( program_string( options.version ) + ' executable file does not exist: ' + 
                  exeFile )
    
    return exeFile, tablesDir

#############################################################################################

def setWorkDirs_and_Flags( inpFiles, silent ):
    """
    Set the paths of the working directories and the inflate flag.
    
    Arguments : inpFiles : list of input files (absolute or relative paths, with or 
                           without extension '.inp' or file name '/INP') 
                silent   : flag indicating that no input must be requested from stdin
       
    Results   : workDirs_and_flags : list of (workDir, inflateFlag) tuples, where
                                     workDir     : working directory
                                     inflateFlag : flag indicating directory has to be 
                                                   inflated; directory does not exist
    
    The input files are either workDir.inp (inflateFlag=True) or workdir/INP (inflateFlag=False)
    """
    
    # determine working directories from input file names
    
    workDirs_and_flags = []
    
    for inpFile in inpFiles:
        
        # absolute path of input file directory
        
        inpDir = os.path.abspath( os.path.dirname( inpFile ) )

        # absolute path of working directory
        
        if os.path.basename( inpFile ) == 'INP':
            workDir = os.path.dirname( inpFile )
        elif os.path.splitext( inpFile )[1] == '.inp':
            workDir = os.path.splitext( inpFile )[0]
        else:
            workDir = inpFile
        
        workDir = os.path.abspath( workDir )
        
        # check whether ".inp" or "/INP" exists and set inflate flag
        
        inpFile1 = workDir + '.inp'
        inpFile2 = os.path.join( workDir, 'INP' )
 
        inflateFlag = ''
        if os.path.isfile( inpFile1 ):
            if os.path.isfile( inpFile2 ):  # both input files exist
                if silent:
                    inflateFlag = 'y'
                else:
                    displayFile( inpFile1 )
                    while inflateFlag != 'y' and inflateFlag != 'n':
                        inflateFlag = input( '\nFile "' + inpFile2 + '" exists.' +
                                             '\nOverwrite file and directory it resides in? ' +
                                             '(y/n): ' )
                        inflateFlag = inflateFlag.lower( )
                if inflateFlag == 'y':
                    print( '\nRemoving directory "' + workDir + '"...' )
                    shutil.rmtree( workDir )
            else:                           # only ".inp" exists
                inflateFlag = 'y'
                if not silent:
                    displayFile( inpFile1 )
        else:
            if os.path.isfile( inpFile2 ):  # only '/INP' exists
                inflateFlag = 'n'
                if not silent:
                    displayFile( inpFile2 )
            else:                           # no input file exists
                sys.exit( '\nNone of the input files "' + inpFile1 + '" and ' +
                          '"' + inpFile2 + '" exists.' )
        
        inflateFlag = (inflateFlag == 'y')
        
        workDirs_and_flags.append( (workDir, inpDir, inflateFlag ) )
        
    return workDirs_and_flags

#############################################################################################

def displayFile( someFile ):
    """
    Display the contents of a file on stdout.
       
    Arguments : someFile - name of the file
       
    Results   : None
    """
    
    f = open( someFile, 'r' )

    print( '\n' + someFile + ':' )
    print( '\n' + f.read( ) )
    
    f.close( )
    
    answer = input( '\nPress <CR> to continue: ' )

#############################################################################################

def getWorkFiles( workDir ):
    """
    Get the work files corresponding to working directory workDir
    """

    workFiles = []
    for workFile in glob.glob( workDir + ".*" ):
        if not os.path.isfile(workFile):          # must be a file to be copied
            continue
        trunk, _ = os.path.splitext( workFile )
        if trunk != workDir:                      # do not copy other fnames containing "."
            continue
        workFiles.append( workFile )
    
    return workFiles

#############################################################################################

def inflate( workDir ):
    """
    Inflate workDir.* into workDir/(*).upper().
    """

    # Create directory workDir unless workDir already exists
    
    if os.path.exists( workDir ):
        print( "WARNING: Path " + workDir + " already exists." )
    else:
        os.mkdir( workDir )
    
    # copy workDir.* to workDir/(*).upper
    
    workFiles = getWorkFiles( workDir )
    for workFile in workFiles:
        _, ext = os.path.splitext( workFile )
        workFile2 = os.path.join( workDir, ext[1:].upper( ) )
        shutil.copy2( workFile, workFile2 )

#############################################################################################

def writePATHfile( workDir, inpDir, tablesDir ):
    """
    Write the file PATH to the directory workDir containing the full path name and the table 
    file path.
    """
    # open PATH file
    
    workDirPath = os.path.join( os.getcwd( ), workDir )        # full directory path name
    pathFile = os.path.join( workDirPath, 'PATH' )
    
    try:
        f_path = open( pathFile, 'w' )
    except:
        sys.exit( 'Error during opening ' + pathFile )

    # write PATH file
    
    f_path.write( '&path' )
    f_path.write( "\nhostname='" + socket.gethostname( ) + "'" )
    f_path.write( "\nworkpath='" + workDirPath + os.sep + "'" )
    f_path.write( "\ninppath='" + inpDir + os.sep + "'" )
    f_path.write( "\ntablepath='" + tablesDir + os.sep + "'" )
    f_path.write( "\n/" )
    
    # close PATH file
    
    f_path.close( )

#############################################################################################

def deflate( workDir, options=None ):
    """
    Move workDir/XXX into workDir.xxx where xxx=XXX.lower().
    
    Only all-uppercase files are moved. If workDir is empty after these moves, it is removed.
    """

    # Check if directory exists
    
    if not os.path.exists( workDir ):
        sys.exit( "Cannot deflate " + workDir + ". Path does not exist." )
    
    # copy workDir/* to workDir.** where **=/(*).upper
    
    deleteDir = True
    for workFile in os.listdir( workDir ):
        
        if os.path.isdir( workFile ) or not workFile.isupper():
            deleteDir = False
            continue
        
        workFile1 = os.path.join( workDir, workFile )
        workFile2 = os.path.join( workDir, os.pardir, 
                              os.path.basename( workDir ) + '.' + workFile.lower( ) )
        workFile2 = os.path.realpath( workFile2 )
        
        if options:
            if options.logfile:
                logfile = os.path.realpath( os.path.abspath( options.logfile ) )
                if workFile2 == logfile:
                    os.remove( workFile1 )          # delete old log file
                    continue

        if workFile == 'PATH':
            os.remove( workFile1 )                  # delete PATH workFile
            continue
        
        shutil.move( workFile1, workFile2 )
    
    if deleteDir:
        os.rmdir( workDir )
    
#############################################################################################

def runProgram( workDir_and_flag, exeFile, tablesDir, options ):
    """
    Run PROGRAM in workDir.
    
    Arguments : exeFile  : absolute path of PROGRAM executable
                workDir  : absolute path of working directory
                options : command line options
       
    Results   : None
    """
    
    # redirect stdout and stderr to logfile if given
    
    if options.logfile:
        logFile = os.path.abspath(options.logfile)
        f_log = open(logFile, 'a')
    else:
        f_log = sys.stdout

    # remember current directory, inflate directory and write PATH file, and change to working 
    # directory
    
    startDir = os.getcwd( )
    
    workDir, inpDir, inflateFlag = workDir_and_flag
    if inflateFlag:
        inflate( workDir )
    writePATHfile( workDir, inpDir, tablesDir )    
    
    os.chdir( workDir )

    # report starting of PROGRAM with date

    date = time.strftime( '%d-%b-%Y %H:%M:%S' )
    f_log.write( program_string( options.version ) + ' started on ' + date + 
                 ' in directory ' + socket.gethostname( ) + ':' + workDir + '\n')

    # apply niceness

    if options.niceness and os.name == 'posix':
        os.nice( options.niceness )
    
    # set the number of threads per job (only working in PROGRAM is compiled with ifort
    # using -qopenmp or gfortran using -fopenmp)
    
    os.environ["OMP_NUM_THREADS"] = str(options.threads)

    # run PROGRAM and postprocessing script POST
    
    f_log.flush( )
    
    try:
        if options.debug:
            if os.name == 'posix':
                subprocess.call( [debugcmd, exeFile], stdout=f_log, stderr=f_log )
            else:
                subprocess.call( [debugcmd, '"' + exeFile + '"'], stdout=f_log, stderr=f_log )
        else:
            subprocess.call( [exeFile], stdout=f_log, stderr=f_log)
        
        date = time.strftime( '%d-%b-%Y %H:%M:%S' )
        f_log.flush( )
        f_log.write( program_string( options.version ) + ' completed on ' + date + 
                     ' in directory ' + socket.gethostname( ) + ':' + workDir + '\n\n')
    except:
        f_log.flush( )
        f_log.write( program_string( options.version ) + 
                     ' terminated abnormally in directory ' + socket.gethostname( ) + ':' + 
                     workDir + '\n\n')
    else:
        if os.path.exists( 'POST' ):
            runpy.run_path( 'POST' )
            date = time.strftime( '%d-%b-%Y %H:%M:%S' )
            f_log.write( '\nPostprocessing completed on ' + date + 
                         ' in directory ' + socket.gethostname( ) + ':' + workDir + '\n')
    
    # change directory back to start directory
    
    os.chdir( startDir )

    # deflate directory (if it was inflated)
    
    if inflateFlag:
        deflate( workDir, options )
        
    # copy files back (if local_host is different from current host)

    local_host = options.localhost
    current_host = socket.gethostname( )
    
    if local_host != current_host:
        if inflateFlag:
            workFiles = getWorkFiles( workDir )
            for workFile in workFiles:
                src = workFile
                dest = local_host + ':' + workFile
                command = 'until rsync -a ' + src + ' ' + dest + '; do sleep 60; done'
                ret = subprocess.call( command, shell=True )
                if ret: 
                    f_log.write( 'Error while copying ' + workFile + ' to ' + local_host + 
                                 '\n' )
        else:
            src = workDir + os.sep
            dest = local_host + ':' + workDir
            command = 'until rsync -a ' + src + ' ' + dest + '; do sleep 60; done'
            ret = subprocess.call( command, shell=True )
            if ret: 
                f_log.write( 'Error while copying "' + workDir + os.sep + '" to ' + local_host + 
                             '\n' )
           
    if f_log is not sys.stdout:
        f_log.close( )
    
#############################################################################################

def main( args ):
    """
    Main function controlling PROGRAM execution.
       
    Arguments : args : command line arguments
       
    Results   : None
    """
    
    # parse the command line
    
    options = parseCommandLine( args, '.inp', '/INP' )
    
    # define various paths
    
    os.chdir( options.startdir )
    exeFile, tablesDir = setProgramPaths( options )    
    workDirs_and_flags = setWorkDirs_and_Flags( options.inpFiles, options.silent )

    # prepare and start script on remote host

    remote_host = options.remotehost
    current_host = socket.gethostname()
    
    if remote_host != current_host:

        # copy tables, executable, and script
        for path in (tablesDir + os.sep, exeFile, args[0]):
            ret = subprocess.call( ['ssh', remote_host, 'mkdir -p ' + os.path.dirname( path )] )
            if ret:
                sys.exit( 'Error while creating ' + os.path.dirname( path ) + ' on ' + 
                          remote_host )
            ret = subprocess.call( ['rsync', '-a', path, remote_host + ':' + path] )
            if ret: 
                sys.exit( 'Error while copying ' + path + ' to ' + remote_host )

        # copy input files
        for workDir, inpDir, inflateFlag in workDirs_and_flags:
            if inflateFlag:
                path = os.path.dirname( workDir )
                ret = subprocess.call( ['ssh', remote_host, 'mkdir -p ' + path] )
                if ret: 
                    sys.exit( 'Error while creating ' + path + ' on ' + remote_host )
                workFiles = getWorkFiles( workDir )
                for workFile in workFiles:
                    ret = subprocess.call( ['rsync', '-a', workFile, 
                                           remote_host + ':' + workFile] )
                    if ret: 
                        sys.exit( 'Error while copying ' + workFile + ' to ' + remote_host )
                    
            else:
                ret = subprocess.call( ['ssh', remote_host, 'mkdir -p ' + workDir] )
                if ret: 
                    sys.exit( 'Error while creating ' + workDir + ' on ' + remote_host )
                ret = subprocess.call( ['rsync', '-a', workDir + os.sep, 
                                        remote_host + ':' + workDir + os.sep] )
                if ret != 0:
                    sys.exit( 'Error while copying "' + workDir + os.sep + '" to ' + 
                              remote_host )

        # start imsil on remote host
        command = '. /opt/intel/oneapi/setvars.sh > /dev/null; '
        if options.background:
            command += 'nohup ' + ' '.join( args ) + ' --local ' + options.localhost + \
                       ' --startdir ' + options.startdir + ' --nobackground -s >> ' + \
                       os.path.join( options.startdir, 'nohup.out' ) + ' 2>&1 &'
            subprocess.Popen( ['ssh', remote_host, command] )
            print( '\n' + program_string( options.version ) + ' has been delegated to ' +
                   remote_host + '\n' )
            sys.exit( )
        else:
            command += ' '.join( args ) + ' --local ' + options.localhost + \
                       ' --startdir ' + options.startdir + ' -s'
            ret = subprocess.call( ['ssh', remote_host, command] )
            if ret: 
                sys.exit( 'Error while submitting ' + program_string( options.version ) + 
                          ' to ' + remote_host )
    else:
    
        # Start script in the background
        
        if options.background:
            command = ['nohup'] + args + ['--nobackground', '-s']
            subprocess.Popen( command )#, stdout=f_null, stderr=f_null ) 
                # cannot use call() here since it waits for the command to complete
                # redirect stdout and stderr to /dev/null to avoid message during wait
            print( '\n' + program_string( options.version ) +
                   ' has been started in the background.\n' )
            sys.exit( )
        
        # run program 
        
        runProgramPartial = functools.partial(runProgram, exeFile=exeFile, tablesDir=tablesDir,
                                              options=options)
        if options.parallel_jobs == 0:
            options.parallel_jobs = len( workDirs_and_flags )
        pool = multiprocessing.Pool( options.parallel_jobs )
        pool.map( runProgramPartial, workDirs_and_flags )
            
    #  display output and surface file
    
    if not options.silent:
        for workDir, inpDir, inflateFlag in workDirs_and_flags:
            if inflateFlag:
                outFile = workDir + '.out'
                surfFile = workDir + '.surf'
            else:
                outFile = os.path.join( workDir, 'OUT' )
                surfFile = os.path.join( workDir, 'SURF' )
            
            if os.path.isfile( outFile ):
                answer = input( '\nDisplay ' + outFile + '? (y/N): ')
                answer = answer.lower( )
                if answer == 'y':
                    displayFile( outFile )
            else:
                print( '\nNo outputfile ' + outFile + ' has been written.' )

            if os.path.isfile( surfFile ):
                answer = input( '\nPlot ' + surfFile + '? (y/N): ')
                answer = answer.lower( )
                if answer == 'y':
                    import plot
                    plot.plot( 'surf', [surfFile] )

#############################################################################################

if __name__ == '__main__':
    main( sys.argv )
#else:
#    print( 'imsil.py loaded as a module' )
