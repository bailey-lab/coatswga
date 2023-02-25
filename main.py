import find as find
import filter as filter
import sets as sets
import numpy as np
import pandas as pd
import sys
import json

def main():
    in_json = sys.argv[2]
    step = sys.argv[1]

    # step = "sets"
    # in_json = '/Users/kaleb/Desktop/Bailey_Lab/code/newswga/new_src/my_params.json'

    with open(in_json, 'r') as f:
        data = json.load(f)
    if step == "find":
        find.main(data)
    elif step == "sets":
        df, primers_with_positions = filter.main(data)
        sets.main(df, primers_with_positions, data)
    else:
        print("Invalid step. Input is 'find' or 'sets'.")

if __name__ == "__main__":
    main()
    
