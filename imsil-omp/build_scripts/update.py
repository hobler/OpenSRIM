import os
import shutil
import subprocess

# These parameters need to be adapted if the script is used for a different
# program:
########################################################################
program = 'imsil'        # name of the program
########################################################################
# This parameter determines where the program is built:
########################################################################
root_dir = os.path.expanduser('~') + '/work/' + program   # root of program workspace
########################################################################
srcdir = os.getcwd()
#print('srcdir: ', srcdir)


#command = 'git describe --tags --exact-match || ' + \
#          'git symbolic-ref -q --short HEAD'
#command = 'if test "`git status -z`"; then ' + \
#          'git symbolic-ref --short HEAD; else ' + \
#          'git describe --tags --exact-match || git symbolic-ref --short HEAD; fi'
#p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
p = subprocess.Popen('git symbolic-ref --short HEAD', 
                     shell=True, stdout=subprocess.PIPE)
branch = p.stdout.read()[:-1].decode()   # strip newline
print('branch:', branch)


def update(project):
    """
    Copy files under the project directory of the git repository to the 
    workspace, and create necessary directories.
    
    Argument:
        project (string): Name of the project (e.g., 'code', 'scripts', etc.)
    """
    print('Running update on ', project)
    proj_dir = os.path.join(root_dir, branch, project)
    if not os.path.isdir(proj_dir):
        os.mkdir(proj_dir)
        #print 'proj_dir: ', proj_dir
        print('Directory ' + proj_dir + ' has been created.')
    
    for dirname, subdirnames, filenames in os.walk(srcdir):
        for filename in filenames:
            file_, ext = os.path.splitext(filename) 
            srcfile = os.path.relpath(os.path.join(dirname, filename))
            
            dstfile = os.path.join(proj_dir, srcfile)
            if dstfile.endswith('~'):
                continue
            
            dstdir = os.path.dirname(dstfile)
            if not os.path.isdir(dstdir):
                os.makedirs(dstdir)
    
            if os.path.exists(dstfile):
                # rounding is necessary since copy2 introduces small errors
                # into the modification time
                src_mtime = round(os.path.getmtime(srcfile), 3)
                dst_mtime = round(os.path.getmtime(dstfile), 3)
                if dst_mtime > src_mtime:
                    print('WARNING: ' + filename + \
                          ' is newer in working directory than in git repository.')
                copy = dst_mtime < src_mtime
            else:
                copy = True
            
            if copy:
                print('Copying ' + srcfile + ' to ' + dstfile)
                shutil.copy2(srcfile, dstfile)

    print('Project has been updated.')
