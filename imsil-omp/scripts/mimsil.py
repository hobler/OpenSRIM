#! /usr/bin/env python3
"""
Usage: python3 mimsil.py [options] [input-file ...]

The purpose of this script is to run imsil, in the following referred to as PROGRAM,
repeatedly (multiple execution). The input files '.inp' are generated from a template 
file '.tem' by replacing place holders '#...#' by values found in a parameter file 
'.par'. The '.par' file contains the names of the input files, and for each input file 
key=value pairs where '#key#' is to be replaced by 'value' in the '.tem' file.

options    : -v <version>       : Version of PROGRAM to be used (name of the subdirectory of 
                                  the PROGRAM code directory). Default is 'trunk'
             -p <parallel_jobs> : Maximum number of jobs to run in parallel
                                  parallel_jobs=0 means run all jobs in parallel
                                  (default is parallel_jobs=1)
             -t <threads>       : Number of threads (coarray images) per job
                                  (default is threads=1)
             -r <hostname>      : host to run PROGRAM on
             -n <niceness>      : nice command for Unix. <niceness> means the niceness value 
                                  (generally between 1 and 20)
             -l <logfile>       : logfile, replacing standard output or nohup.out
             -s                 : silent execution (no questions asked)
             -i                 : interactive (foreground) execution
             -d                 : debug
             -D                 : use code with debug information (default if -d is set)
             -O                 : use optimized code (default if -d is not set)

             If -i is not set and -l is not set, the logfile is set as the stem of the first
             .tem file + '.log'. If logfile exists, it is overwritten.
             
             Unless -s is set, the user is prompted for modifying these options.

input-file : the name of the template or parameter file. The extension '.tem' or '.par'
             may be appended or omitted.
"""

import collections, os, io, sys, socket, time
import imsil

#############################################################################################
def program_string( version ):
    """
    Return program and version information as a string.
    Meant for use in output messages.
    """
    return 'MIMSIL(' + version + ')'

#############################################################################################

def setPaths( temFile, options ):
    """
    Make input files from template and parameter file.
       
    Argument : temFile : template or parameter file (with or without extension,
                         absolute or relative path)
               options : Options
    
    Results  : temFile : absolute path of template file
               parFile : absolute path of parameter file
    """
    
    # strip extension and get template and parameter file name
    
    if os.path.splitext( temFile )[1] in ['.tem', '.par']:
        trunk = os.path.splitext( temFile )[0]
    else:
        trunk = temFile
        
    trunk = os.path.abspath( trunk )
    
    temFile = trunk + '.tem'
    parFile = trunk + '.par'
    ptemFile = trunk + '.ptem'
    
    # display parameter and template file

    if os.path.isfile( parFile ):
        if not options.silent:
            imsil.displayFile( parFile )
    else:
        sys.exit( 'Parameter file ' + parFile + ' does not exist' )
    
    if os.path.isfile( temFile ):
        if not options.silent:
            imsil.displayFile( temFile )
    else:
        sys.exit( 'Template file ' + temFile + ' does not exist' )
    
    if os.path.isfile( ptemFile ):
        if not options.silent:
            imsil.displayFile( ptemFile )

    return( temFile, parFile )

#############################################################################################

def makeInpFiles( temFile ):
    """
    Make input files from template and parameter file.
       
    Arguments : temFile : absolute path of template or parameter file 
                          (extension .tem or .par optional)
       
    Results  : List of input files
    """
    
    # read template and parameter file
    
    fname, ext = os.path.splitext( temFile )
    temFile = fname + '.tem'
    parFile = fname + '.par'
    ptemFile = fname + '.ptem'
    
    f_tem = open( temFile, 'r' )
    temString = f_tem.read( )
    f_tem.close( )
    
    f_par = open( parFile, 'r' )
    parString = f_par.read( ).strip( )
    f_par.close( )
    
    if os.path.exists( ptemFile ):
        f_ptem = open( ptemFile, 'r' )
        ptemString = f_ptem.read( )
        f_ptem.close( )
    else:
        ptemString = None 
    
    rootDir = os.path.dirname( temFile )

    # parse the parameter file and define the param dictionary
    
    param = collections.OrderedDict( )
    
    if parString.startswith( 'FILE' ):
        
        for line in io.StringIO( parString ):
            if not line.strip( ):                # skip empty lines
                continue
            elif line.strip( ).lower( ) == 'end': # premature end
                break
            elif line.startswith( 'FILE:' ):
                inpFile = line[5:].strip( )
                param[inpFile] = {}
            else:
                key_and_value = line.split( '=', 1 )
                if len( key_and_value ) < 2:
                    sys.exit( 'Error reading the following line from the parameter file:'
                              '\n' + line )
                key, value = key_and_value
                key = key.strip( )
                value = value.strip( )
                param[inpFile][key] = value

    elif parString.startswith( 'def ' ):
        
        from types import FunctionType
        fun_code = compile( parString, '<string>', 'exec' )
        get_params = FunctionType( fun_code.co_consts[0], globals(), 'get_params')
        get_params( param )
    
    else:
        
        exec( parString )

    # substitute values for keys
    
    inpFiles = list(param.keys( ))
    
    for inpFile in inpFiles:
        inpString = temString
        for key, value in list(param[inpFile].items( )):
            inpString = inpString.replace( '#' + key + '#', value )
        f_inp = open( os.path.join( rootDir, inpFile + '.inp' ), 'w' )
        f_inp.write( inpString )
        f_inp.close( )
    
    if ptemString is not None:
        for inpFile in inpFiles:
            postString = ptemString
            for key, value in list(param[inpFile].items( )):
                postString = postString.replace( '#' + key + '#', value )
            f_post = open( os.path.join( rootDir, inpFile + '.post' ), 'w' )
            f_post.write( postString )
            f_post.close( )
    
    inpFiles = [os.path.join( rootDir, inpFile + '.inp') for inpFile in inpFiles]
    
    return inpFiles

#############################################################################################

def makeInpFilesOld( temFile ):
    """
    Make input files from template and parameter file.
       
    Arguments : temFile : absolute path of template or parameter file 
                          (extension .tem or .par optional)
       
    Results  : List of input files
    """
    
    # read template file
    
    fname, ext = os.path.splitext( temFile )
    temFile = fname + '.tem'
    parFile = fname + '.par'
     
    f_tem = open( temFile, 'r' )
    temString = f_tem.read( )
    f_tem.close( )
    
    rootDir = os.path.dirname( temFile )

    # parse the parameter file and substitute values for keys
    
    inpFile = ''
    inpFiles = []
    f_par = open( parFile, 'r' )

    for line in f_par:
        if not line.strip( ):                # skip empty lines
            continue
        elif line.strip( ).lower( ) == 'end': # premature end
            break
        elif line.startswith( 'FILE:' ):
            if inpFile:                     # write previous input file to disk
                # check for "#.*" in input
                f_inp = open( inpFile, 'w' )
                f_inp.write( inpString )
                f_inp.close( )
            inpFile = os.path.join( rootDir, line[5:].strip( ) ) + '.inp'
            inpFiles.append( inpFile )
            inpString = temString
        else:
            key_and_value = line.split( '=', 1 )
            if len( key_and_value ) < 2:
                sys.exit( 'Error reading the following line from the parameter file:'
                          '\n' + line )
            key, value = key_and_value
            key = key.strip( )
            value = value.strip( )
            inpString = inpString.replace( '#' + key + '#', value )
                                            
    if inpFile:                             # write last input file to disk
        f_inp = open( inpFile, 'w' )
        f_inp.write( inpString )
        f_inp.close( )
    
    f_par.close( )
    
    return inpFiles

#############################################################################################


def main( ):
    """
    Main function controlling the multiple execution of PROGRAM.
    """
    
    # parse the command line
    
    options = imsil.parseCommandLine( sys.argv, '.tem', '.par' )
    temFiles = options.inpFiles
    
    # set the paths of template and parameter file and generate input files
    
    inpFiles = []
    
    for temFile in temFiles:
        temFile, parFile = setPaths( temFile, options )
        inpFiles = inpFiles + makeInpFiles( temFile )
    
    # call imsil with input files
    
    programpath = os.path.join( os.path.dirname( os.path.realpath( sys.argv[0] ) ), 'imsil.py' ) 
    args = [ programpath, '-v', options.version, '-n', str(options.niceness), 
            '-p', str(options.parallel_jobs), '-t', str(options.threads), 
            '-r', options.remotehost, '-s' ] 
    if not options.background:
        args.append( '-i' )
    if options.debug:
        args.append( '-d' )
    if options.debugcode:
        args.append( '-D' )
    if options.logfile:
        args.append( '-l')
        args.append( options.logfile )
#    if (options.background or options.silent) and not options.logfile:
#        args.append( '-l')
#        logfile = os.path.splitext( temFiles[0] )[0] + '.log'
#        args.append( logfile )
#    else:
#        logfile = options.logfile
    args = args + inpFiles

#    if options.logfile or options.background or options.silent:
    if options.logfile:
        if os.path.isfile( options.logfile ):
            os.remove( options.logfile )
    
    date = time.strftime( '%d-%b-%Y %H:%M:%S' )
    print( program_string( options.version ) + ' started on ' + date + 
           ' with files ', temFiles)

    imsil.main( args )

    if not options.background:
        date = time.strftime( '%d-%b-%Y %H:%M:%S' )
        print( program_string( options.version ) + ' completed on ' + date + 
               ' with files ', temFiles, '\n')

#############################################################################################

if __name__ == '__main__':
    main( )
else:
    print( 'mimsil.py loaded as a module' )