import json
import sys
import find as find
import filter as filter
import sets as sets
from datetime import datetime

def main():
    if len(sys.argv) != 3:
        print('Usage: python3 main.py <STEP: find|sets|all> [PARAMS FILEPATH]')
        quit()
    else:
        in_json = sys.argv[2]
        step = sys.argv[1]

        with open(in_json, 'r') as f:
            data = json.load(f)
        
        if not data['data_dir'][-1] == "/":
            data['data_dir'] = data['data_dir'] + "/"

        print("=" * 80)
        for elt in data:
            print(f'{elt}: {data[elt]}')
        dt = datetime.now().strftime('%B %d, %Y at %H:%M')
        print(f"Started on {dt}")
        print("-" * 80)

        if step == "find":
            find.main(data)
        elif step == "sets":
            df, primers_with_positions, chr_lens = filter.main(data)
            sets.main(df, primers_with_positions, chr_lens, data)
        elif step == "all":
            find.main(data)
            df, primers_with_positions, chr_lens = filter.main(data)
            sets.main(df, primers_with_positions, chr_lens, data)
        else:
            print("Invalid step. Input is 'find', 'sets', or 'all'.")

if __name__ == "__main__":
    # in_json = '/Users/kaleb/Desktop/Bailey_Lab/code/newswga/params/new_params.json'
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
    
