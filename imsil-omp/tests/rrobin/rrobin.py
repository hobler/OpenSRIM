print('__file__:', __file__)

if __file__ == 'plot.py':
    from rrobin.rrobin_fun import rrobin_fun
else:
    from rrobin_fun import rrobin_fun

rrobin_fun('rrobin')
