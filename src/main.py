import json
import sys
from . import find
from . import filter
from . import sets
from datetime import datetime
import argparse

defaults = {
    "write": False,
    "verbose": False,
    "cpus": 1,
    "min_primer_length": 8,
    "max_primer_length": 16,
    "min_tm": 15,
    "max_tm": 45,
    "min_fg_count": 200,
    "max_ratio": 1,
    "fragment_length": 2000,
    "target_coverage": 0.95,
    "force_coverage_threshold": False,
    "min_set_size": 10,
    "existing_primers": [],
}

class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


def main():
    fmt = lambda prog: CustomHelpFormatter(prog)
    par = argparse.ArgumentParser(prog='swga3',description="Finds selective primers that preferentially bind to the target genome.", formatter_class=fmt)
    par.add_argument('-j', '--json_file', metavar='<path>', help='Path to the file containing parameters in JSON format')
    par.add_argument('-d', '--data_dir', metavar='<path>', help='Path to the directory to store kmer and CSV files')
    par.add_argument('-v', '--verbose', action='store_true', help='Whether or not extra output is wanted (intermediate sets, timing of parts), default: False')
    par.add_argument('-w', '--write', action='store_true', help='Whether or not CSV\'s of the primers with counts should be written')
    par.add_argument('-bg', '--bg_genomes', action='extend', nargs='*', metavar='<path1> <path2>', help='file paths of the background (off-target) genomes')
    par.add_argument('-fg', '--fg_genomes', action='extend', nargs='*', metavar='<path1> <path2>', help='file paths of the foreground genomes')
    par.add_argument('-bp', '--bg_prefix', help='Prefix for the background files')
    par.add_argument('-fp', '--fg_prefixes', action='extend', nargs='*', help='List of prefixes for the foreground files')
    par.add_argument('-c', '--cpus', type=int, help="Number of CPUs to use")
    par.add_argument('-m', '--min_primer_length', type=int, help="Minimum primer length")
    par.add_argument('-M', '--max_primer_length', type=int, help="Maximum primer length")
    par.add_argument('-t', '--min_tm', type=int, help="Minimum predicted primer melting temperature")
    par.add_argument('-T', '--max_tm', type=int, help="Maximum predicted primer melting temperature")
    par.add_argument('-g', '--min_fg_count', type=int, help="Minimum number of foreground genome occurrences, helps to filter out less effective primers")
    par.add_argument('-r', '--max_ratio', type=int, help="Maximum value of the background occurrences divided by foreground occurrences")
    par.add_argument('-l', '--fragment_length', type=int, help="Theoretical length of DNA fragments formed by DNA polymerase enzymes")
    par.add_argument('-o', '--target_coverage', type=int, help="Proportion of the genome to cover")
    par.add_argument('-f', '--force_coverage_threshold', action='store_true', help="Force the program to run until either every potential primer has been checked or the target coverage threshold has been reached")
    par.add_argument('-s', '--min_set_size', type=int, help="Minimum number of primers to include in the set")
    par.add_argument('-p', '--existing_primers', action='extend', nargs='*', metavar='<primer1> <primer2>', help='set of primers to check against, defualt: None')
    args = par.parse_args()
    args = vars(args)

    if args['json_file']:
        with open(args['json_file'], 'r') as f:
            data = json.load(f)
        if data['write'] and not args['write']:
            args['write'] = True
        if data['verbose'] and not args['verbose']:
            args['verbose'] = True
        if data['force_coverage_threshold'] and not args['force_coverage_threshold']:
            args['force_coverage_threshold'] = True
    else:
        data = {}
    del args['json_file']
    
    breaker = False
    for arg in args:
        if args[arg] is not None:
            data[arg] = args[arg]
        if arg not in data:
            if arg in defaults:
                data[arg] = defaults[arg]
            else:
                print(f"{arg} not specified")
                breaker = True
    if breaker:
        quit()

    if data['data_dir'][-1] != "/":
        data['data_dir'] = data['data_dir'] + "/"

    print("=" * 80)
    for elt in data:
        print(f'{elt}: {data[elt]}')
    dt = datetime.now().strftime('%B %d, %Y at %H:%M')
    print(f"Started on {dt}")
    print("-" * 80)

    find.main(data)
    df, primers_with_positions, chr_lens = filter.main(data)
    sets.main(df, primers_with_positions, chr_lens, data)

if __name__ == "__main__":
    # in_json = '/Users/kaleb/Desktop/Bailey_Lab/code/newswga/params/test_params.json'
    # with open(in_json, 'r') as f:
    #     data = json.load(f)
    
    # if not data['data_dir'][-1] == "/":
    #     data['data_dir'] = data['data_dir'] + "/"
    
    # print("=" * 80)
    # for elt in data:
    #     print(f'{elt}: {data[elt]}')
    # print("-" * 80)
    
    # find.main(data)
    # df, primers_with_positions, chr_lens = filter.main(data)
    # sets.main(df, primers_with_positions, chr_lens, data)
    main()
    
