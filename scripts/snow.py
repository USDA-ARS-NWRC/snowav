
import sys
import snowav

'''

'''

def run(config_file):

    snow = snowav.framework.framework.SNOWAV(config_file = config_file)

if __name__ == '__main__':
    config_file = sys.argv[1]
    run(config_file)
