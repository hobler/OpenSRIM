#!/usr/bin/python3
"""
Search for dependencies in IMSIL and write files for inclusion in a makefile.

IMSIL is assumed to be under git control. This script has to be called 
from somewhere within the git workspace.

The intended directory structure is as follows:

    src_dir                           .. defined below in this script
        > casc                        \
        > eval                         > containing fortran files defining
        > ...                         /     modules
        > unittests
            > casc                    \
            > eval                     > containing unittest main progams and 
            > ...                     /     *.mk files for make of unittests
        main.f90                      .. main program of imsil
        makefile
        *.mk                          .. include files for make of imsil.exe

    exe_dir                           .. defined below in this script
        > obj                         .. containing .o files
        > obj_db
        > mod                         .. containing .mod files
        > mod_db
        imsil.exe                     .. IMSIL executable
        imsil_db.exe
        > unittests
            > casc                    \
            > eval                     > containing test executables
            > ...                     / 

Output Files:
    dirs.mk, dirs_db.mk: Macros SRCDIR, OBJDIR, MODDIR, EXE, defining the root 
        directories for the source and object files, the directory for mod 
        files, and the path of the executable, respectively. In addition,
        makes all necessary directories
    objs.mk: Macro OBJS, defining the paths of all object files relative to
        OBJDIR.
    depend.mk: Dependency lines and actions to be performed.
    
The output files are only written if they differ from the previous version. This
avoids makefile to remake the target if nothing has changed.
"""
#from future.backports.test.pystone import TRUE
from numpy.f2py.crackfortran import endifpattern
from _curses import version

USAGE = 'Usage: python make_code.py main_relpath make_args'

import sys
import os
import filecmp
import re
import subprocess

# Root directories and program names for source and executable code
########################################################################
# root directory of the source code
src_dir = os.path.join(os.path.expanduser('~'), 'git/imsil/code')
# directory for executable code and all auxiliary files (.o, .mod, make 
# include files)
p = subprocess.Popen('git symbolic-ref --short HEAD', 
                     shell=True, stdout=subprocess.PIPE)
branch = p.stdout.read()[:-1].decode("utf-8")   # strip newline
#print('Type of branch:', type(branch))
print('branch:', branch)
exe_dir = os.path.join(os.path.expanduser('~'), 'work/imsil', branch, 'code')
# name of the file (without extension) containing the main  program source 
# relative to src_dir
main_relpath = sys.argv[1] + '.f90'
main_reldir = os.path.dirname(main_relpath)
main_src_dir = os.path.join(src_dir, main_reldir)
main_exe_dir = os.path.join(exe_dir, main_reldir)
# Fortran predefined modules (to avoid warnings)
predefined_modules = ['f95_precision', 'lapack95', 'iso_fortran_env']
########################################################################
# Here is some untested code to amend the directory currently named after 
# the git branch by the tag of the revision used to create the code.
#command = 'git describe --tags --exact-match || ' + \
#          'git symbolic-ref -q --short HEAD'
#command = 'if test "`git status -z`"; then ' + \
#      'git symbolic-ref --short HEAD; else ' + \
#      'git describe --tags --exact-match || git symbolic-ref --short HEAD; fi'
#p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

nwarn = 0

def move_file(src, dst):
    """
    Intelligently move file src to dst.
    
    Moves src to dst, if dst does not exist or is different from src.
    Deletes src, if src and dst do not differ.
    
    Arguments:
        src: path of file to be moved
        dst: path to be moved to
    """
    if not os.path.exists(dst):
        os.rename(src, dst)
    else:
        if filecmp.cmp(src, dst):
            os.remove(src)
        else:
            os.remove(dst)
            os.rename(src, dst)


def make_depend_objs():
    """
    Write the make include files depend.mk and objs.mk, containing dependency 
    lines and a macro defining all object files, respectively. 
    
    Determines all files needed by main_relpath by recursively traversing
    modules appearing in use statements. Only considers subdirectories of
    src_dir. Modules not found are ignored, i.e. it is assumed they will be 
    linked in the makefile.
    """

    def get_modules(file_relpath):
        """
        Get the names of the modules defined in the file file_relpath.
        
        Argument:
            file_relpath: path relative to src_dir of the source file to be 
                scanned for module definitions
        
        Returns:
            modules (list): list of modules defined in file_relpath
        """
        modules = []
        with open(os.path.join(src_dir, file_relpath), 'r') as f:
            regexp = re.compile(r'^[ \t]*module[ ]+(?P<name>[a-zA-Z][_0-9a-zA-Z]*)[ \t]*$')
            for line in f.readlines():
                result = regexp.search(line)
                if result == None:
                    continue
                #print 'Found module ' + result.group('name').lower() 
                module = result.group('name').lower()
                modules.append(module)
        
        return modules

    def make_file_dict():
        """
        Construct the mapping between source files and module names. 
        
        The dictionary file_dict is a nonlocal variable defined in the 
        enclosing function. Its keys are the module names, and its values are 
        the names of the files defining the respective module.
         
        All source files found in src_dir and its subdirectories are scanned. 
        """
        global nwarn
        
        for dirname, subdirnames, filenames in os.walk(src_dir):
            for filename in filenames:
                _, ext = os.path.splitext(filename) 
                if ext != '.f90':
                    continue
                file_relpath = os.path.relpath(os.path.join(dirname, filename), 
                                               src_dir)
                
                modules = get_modules(file_relpath)
                
                for module in modules:
                    if module in file_dict:
                        print('WARNING: multiple definition of module ', module)
                        nwarn += 1
                        continue
                    file_dict[module] = file_relpath

    def get_prerequisites(file_relpath):
        """
        Determine a list of prerequisites for file_relpath.
        
        A prerequisite for file_relpath is the file defining a module 
        referenced in a use statement of file_relpath.
        
        Argument:
            file_relpath: relative path of the file to be scanned
        
        Returns:
            prerequisites (list): list of relative paths of prerequisite files
        """
        global nwarn
        prerequisites = []
        
        # open FORTRAN file
        path = os.path.join(src_dir, file_relpath)
        if not os.path.exists(path):
            sys.exit('Missing source file:' + path)
        with open(path, 'r') as f:
            # regular expression to find a use statement
            re_use = re.compile(r'^[ \t]*use[ ]+(?P<name>[a-zA-Z][_0-9a-zA-Z]*).*')
            re_ifdef = re.compile(r'^#ifdef[ ]+(?P<var>[a-zA-Z][_0-9a-zA-Z]*).*')
            # search for use statements and extract the prerequisites
            valid = True
            for line in f.readlines():
                result = re_ifdef.search(line)
                if result is not None:
                    valid = (result.group('var') == VERSION)
                    continue
                if line.strip() == '#else':
                    valid = not valid
                    continue
                if line.strip() == '#endif':
                    valid = True
                    continue
                if valid:
                    result = re_use.search(line)
                    if result is not None:
                        module = result.group('name').lower()
                        if module not in file_dict:
                            if module not in predefined_modules:
                                print('WARNING: undefined module ' + module)
                                nwarn += 1
                            continue
                        prerequisites.append(file_dict[module])

        return prerequisites

    def get_dependency_line(target, prerequisites):
        """
        Assemble makefile dependency line for the target target.
        
        Arguments:
            target: relative path of target
            prerequisites (list): list of relative paths of prerequisites 
        
        Returns:
            dependency_line: makefile dependency line
        """
        objfile = os.path.join('${OBJDIR}', os.path.splitext(target)[0]+'.o')
        srcfile = os.path.join('${SRCDIR}', target)
        dependency_line = objfile + ' : ' + srcfile
        for prerequisite in prerequisites:
            dependency_line += ' \\\n\t' + os.path.join('${SRCDIR}', 
                                                        prerequisite)
        dependency_line += ' \n\t' + '${FC} ${FFLAGS} -o $@ -c ' + srcfile

        return dependency_line

    def treat_target(target):
        """
        Recursive procedure treating a target:
        
        Add dependency line to dependency line list.
        
        Add target to objects list if no prerequisites exist that are
        not yet in the objects list or if prerequisite equals target (i.e.,
        if module is defined in the same file as it is used; otherwise 
        treat prerequisites that are not in the objects list as targets.
        
        target still has the extension '.f90'. Must be changed to '.o' later.
        
        Argument:
            target: relative file name of target
        """
        #print 'Looking for dependencies in ', target, '...'
        # determine prerequisites
        prerequisites = get_prerequisites(target)
        #print 'prerequisites=', prerequisites
        #print 'objects=', objects
        #raw_input('Enter <cr> to continue')
        try:
            # treat as targets all prerequisites not yet in the object list 
            for prerequisite in prerequisites:
                if prerequisite not in objects and prerequisite != target:
                    treat_target(prerequisite)
        except RuntimeError:
            sys.exit( 'RuntimeError when processing prerequisite ' + prerequisite +
                      '. Circular dependence?')
        
        # add target to object list and dependency line to list
        object = target                      
        if object not in objects:
            objects.append(object)

        dependency_line = get_dependency_line(target, prerequisites)
        if dependency_line:
            dependency_lines.append(dependency_line)
            
    # dictionary defining the mapping between module and file names
    file_dict = {}
    # list of dependency lines to be written to the makefile
    dependency_lines = []
    # list of objects ordered such that list elements depend on previous list elements only
    objects = []                             

    # set up dictionaries
    make_file_dict()

    # produce the object_line and dependency_lines
    print('Looking for dependencies ...', end=' ')
    treat_target(main_relpath)
    dependency_lines = '\n\n'.join(dependency_lines)

    # write file
    tmp_file = os.path.join(main_src_dir, 'depend_tmp.mk')
    perm_file = os.path.join(main_src_dir, 'depend.mk')

    with open(tmp_file, 'w') as f:
        f.write(dependency_lines)
        
    move_file(tmp_file, perm_file)

    if nwarn > 0:
        print('make_depend: ', nwarn, ' WARNINGS')
        #raw_input('\nPress return to exit... ')
    
    make_srcs_list(objects)
    
    make_objs_list(objects)
    
    print('done')

def make_srcs_list(objects):
    """
    Write the make include file srcs.mk, containing the definition of the 
    SRCS macro, for inclusion in a makefile. SRCS defines all the source files.

    Argument:
        objects (list): File names of the source files relative to src
            root directory.
    """
    object_line = 'SRCS = '
    for object in objects:
        object_line += '\\\n\t' + \
                os.path.join('${SRCDIR}', object)

    tmp_file = os.path.join(main_src_dir, 'srcs_tmp.mk')
    perm_file = os.path.join(main_src_dir, 'srcs.mk')

    with open(tmp_file, 'w') as f:
        f.write(object_line)
    
    move_file(tmp_file, perm_file)
    

def make_objs_list(objects):
    """
    Write the make include file objs.mk, containing the definition of the 
    OBJS macro, for inclusion in a makefile. OBJS defines all the object files.

    Argument:
        objects (list): File names of the source files relative to src
            root directory.
    """
    object_line = 'OBJS = '
    for object in objects:
        obj_relpath = os.path.splitext(object)[0]
        object_line += '\\\n\t' + \
                os.path.join('${OBJDIR}', obj_relpath+'.o')

    tmp_file = os.path.join(main_src_dir, 'objs_tmp.mk')
    perm_file = os.path.join(main_src_dir, 'objs.mk')

    with open(tmp_file, 'w') as f:
        f.write(object_line)
    
    move_file(tmp_file, perm_file)
    

def make_macros(db, omp):
    """
    Write the file dirs.mk or dirs_db.mk, defining the macros OBJDIR, MODDIR, 
    EXE.
    
    OBJDIR is the root directory of the object file directories. 
    MODDIR is the directory for the module files.
    EXE is the path of the executable. 
    
    Arguments:
        db (bool): Flag indicating compilation with debug options.
        omp (bool): Flag indicating compilation with OpenMP
    """
    obj_dir = 'obj'
    mod_dir = 'mod'
    lib_dir = 'lib'
    
    if main_reldir:
        exe = os.path.splitext(main_relpath)[0]
    else:
        exe = 'imsil'
        
    if db:
        obj_dir += '_db'
        mod_dir += '_db'
        lib_dir += '_db'
        exe += '_db'
    if omp:
        obj_dir += '_omp'
        mod_dir += '_omp'
        lib_dir += '_omp'
        exe += '_omp'
    exe += '.exe'
    
    macro_lines = 'SRCDIR=' + src_dir.replace(os.path.expanduser('~'), '~', 1)
    macro_lines += '\nEXEDIR=' + exe_dir.replace(os.path.expanduser('~'), '~', 1)
    macro_lines += '\n'
    macro_lines += '\nOBJDIR=${EXEDIR}' + os.sep + obj_dir
    macro_lines += '\nMODDIR=${EXEDIR}' + os.sep + mod_dir
    macro_lines += '\nLIBDIR=${EXEDIR}' + os.sep + lib_dir
    macro_lines += '\nEXE=${EXEDIR}' + os.sep + exe

    # make sure all the directories exist including subdirectories
    macro_lines += '\n'
    macro_lines += '\n$(shell mkdir -p ${EXEDIR})'
    macro_lines += '\n$(shell mkdir -p ${OBJDIR})'
    macro_lines += '\n$(shell mkdir -p ${MODDIR})'
    if main_reldir:
        macro_lines += '\n$(shell mkdir -p ${EXEDIR}/' + main_reldir + ')' 
    
    sub_dirs = [os.path.relpath(x[0], src_dir) for x in os.walk(src_dir)]
    sub_dirs.remove('.')

    macro_lines += '\n'
    for sub_dir in sub_dirs:
        macro_lines += '\n$(shell mkdir -p ${OBJDIR}' + os.sep + sub_dir + ')'
    
    tmp_file = os.path.join(main_src_dir, 'dirs_tmp.mk')
    perm_file = 'dirs'
    if db:
        perm_file += '_db'
    if omp:
        perm_file += '_omp'
    perm_file = os.path.join(main_src_dir, perm_file + '.mk')
    
    with open(tmp_file, 'w') as f:
        f.write(macro_lines)
    
    move_file(tmp_file, perm_file)
    

def make_libs_file(db, omp):
    """
    Write the make include file libs.mk for inclusion in a makefile. It defines the macro LIBS 
    as a list of all libraries to be linked.

    Arguments:
        db (bool): Flag indicating compilation with debug options.
        omp (bool): Flag indicating compilation with OpenMP
    """
    library_line = 'LIBS = '
    if VERSION == 'AMAT':
        library_line += '\\\n\t' + \
            os.path.join('${LIBDIR}', 'lib_relax.a')
    elif VERSION == 'GTS':
        library_line += '\\\n\t' + \
            os.path.join('${LIBDIR}', 'libimsil_lshyperoctree.a')
            #os.path.join('${LIBDIR}', 'libimsil_lshyperoctree.so')
    
    with open(os.path.join(main_src_dir, 'libs.mk'), 'w') as f:
        f.write(library_line)
    

def main():
    """
    Main program.
    """
    global VERSION
    VERSION = ''

    # Collect make arguments
    if len(sys.argv) == 2:
        make_args = ''
    else:
        args = []
        for arg in sys.argv[2:-1]:
            if '=' in arg:
                args.append(arg)
                if arg.split('=')[0] == 'VERSION':
                    VERSION = arg.split('=')[1]
                    print('Version:', VERSION)
            else:
                args.append(arg+'=1')   # pass preprocessor variable as make parameter
        args.append(sys.argv[-1])
        make_args = ' '.join(args)
            
    # Create the include files for makefile
    make_depend_objs()
    
    print(80*'*')

    for db, omp in (('db', ''), ('db', 'omp'), ('', ''), ('', 'omp')):
        # Create include file defining macros
        make_macros(db, omp)
        # Create include file defining libraries
        make_libs_file(db, omp)
    
        # Create shell command to start make
        make_cmd = 'make '
        if main_reldir:
            make_cmd += '-f ' + os.path.join(src_dir, 'makefile') + ' '
        if db:
            make_cmd += 'DB=db '
        if omp:
            make_cmd += 'OMP=omp '
        make_cmd += make_args
        print('\nRunning "' + make_cmd + '" ...')
        sys.stdout.flush()
        command = 'cd ' + main_src_dir + '; ' + make_cmd
        #print 'command: ', command
        #print 'PATH: ', os.environ['PATH']

        # Run make and collect output
        p = subprocess.Popen(command, shell=True, #executable='/bin/bash', 
                             #env={'PATH' : '.:'+os.environ['PATH']}, 
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        out = out.decode("utf-8")
        err = err.decode("utf-8")
        
        print("\n" + out)
        if err:
            print("Errors and warnings:\n")
            print(err)
            if 'error #' in err:
                print("Aborting because of severe errors.")
                break
            else:
                print("No severe errors.")
        else:
            print("No errors.")
            
        print(80*'*')

if __name__ == '__main__':
    main()