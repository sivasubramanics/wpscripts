#!/usr/bin/env python3

# This script takes a JSON file and formats it with indentation.

import json

def format_json(input_file, output_file):
    with open(input_file, 'r') as f:
        json_string = f.read().strip()
    parsed_json = json.loads(json_string)
    formatted_json = json.dumps(parsed_json, indent=2)
    with open(output_file, 'w') as f:
        f.write(formatted_json)

# Example usage
# Usage: python3 json_formatter.py input.json output.json
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print('Usage: python3 json_formatter.py input.json output.json')
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    format_json(input_file, output_file)
