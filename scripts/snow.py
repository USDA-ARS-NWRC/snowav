
import sys
import snowav

'''
Issues:
- make sure there are unique records; when 'overwrite', should be deleting
current records and then replacing
- memory issues?
- run metadata table
- store elevation bands in array

'''

def run(config_file):

    snow = snowav.framework.framework.SNOWAV(config_file = config_file)

if __name__ == '__main__':
    config_file = sys.argv[1]
    run(config_file)
