import os, sys

if len(sys.argv) == 2:
    old_path = os.path.join('..', '..', sys.argv[1], 'tests')
    if not os.path.exists(old_path):
        sys.exit(f'old path {old_path} does not exist.')
else:
    old_path = os.path.join('..', 'tests_ref')
    if not os.path.exists(old_path):
        old_path = ''

folders = sorted(os.listdir('.'))
for folder in folders:
    if not os.path.isdir(folder):
        continue
    out_files = sorted(os.listdir(folder))
    for out_file in out_files:
        if not out_file.endswith('.out') or out_file in ('nohup.out', 'gmon.out'):
            continue
        
        out_path = os.path.join(folder, out_file)
        last_line = ''
        with open(out_path) as f:
            for line in f:
                last_but_one_line = last_line
                last_line = line
        if not last_but_one_line:
            continue
        try:
            exetime = last_but_one_line.split()[2]
            walltime = last_line.split()[3]
        except IndexError:
            print(f'{out_file:21}: ERROR in file')
            continue
        
        h, m, s = exetime.split(':')
        if h == '00':
            exetime = f'{int(m):5}:{s}'
        else:
            exetime = f'{int(h):2}:{m}:{s}'
        cexetime = (int(h)*60 + int(m))*60 + int(s)

        h, m, s = walltime.split(':')
        if h == '00':
            walltime = f'{int(m):5}:{s}'
        else:
            walltime = f'{int(h):2}:{m}:{s}'

        # corresponding file in the old_path
        if old_path:
            out_path = os.path.join(old_path, folder, out_file)
            if not os.path.exists(out_path):
                continue
            last_line = ''
            with open(out_path) as f:
                for line in f:
                    last_but_one_line = last_line
                    last_line = line
            old_exetime = last_but_one_line.split()[2]
            
            h, m, s = old_exetime.split(':')
            if h == '00':
                old_exetime = f'{int(m):5}:{s}'
            oexetime = (int(h)*60 + int(m))*60 + int(s)
            if oexetime != 0:
                dexetime = cexetime / oexetime - 1
            else:
                dexetime = 0.
                    
        if not old_path:
            print(f'{out_file:21}: {exetime} {walltime}')
        else:
            if len(sys.argv) == 2:
                print(f'{out_file:21} current: {exetime}  '
                      f'{sys.argv[1]}: {old_exetime}  ({dexetime:+.2%})')
            else:
                print(f'{out_file:21} current: {exetime}  '
                      f'tests_ref: {old_exetime}  ({dexetime:+.2%})')
