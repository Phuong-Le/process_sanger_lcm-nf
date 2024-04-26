import sys
import os
import optparse
from pathlib import Path
from collections import defaultdict

def read_samples_list(file_path):
    sample_dict = {}
    try:
        with open(file_path, 'r') as file:
            # Read the header line and split it into column names
            headers = file.readline().strip().split('\t')
            for line in file:
                # Split each line into values using tab as delimiter
                values = line.strip().split('\t')
                # Create a dictionary where keys are column headers and values are corresponding values from the line
                entry = dict(zip(headers, values))
                # Assuming 'sample_id' is unique and using it as the key in the dictionary
                sample_dict[entry['sample_id']] = entry
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return sample_dict


def read_contamination_file(file_path):
    contamination_dict = {}
    try:
        with open("/lustre/scratch126/casm/team267ms/al35/hackathon2024/output/work/tmp/59/5fed489c809e7f612e79449a5fac4d/contamination.txt", 'r') as file:
            for line in file:
                line = line.strip()
                if line:
                    # Split each line into values using whitespace as delimiter
                    values = line.strip().split('\t')
                    # Assuming 'sample_id' and 'normal_id' are the first two values in the line
                    sample_id = values[0]
                    normal_id = values[1]
                    # Assuming 'sample_contamination' and 'normal_contamination' are the third and fourth values respectively
                    sample_contamination = float(values[2])
                    normal_contamination = float(values[3])
                    contamination_dict[sample_id] = {'contamination': sample_contamination}
                    contamination_dict[normal_id] = {'contamination': normal_contamination}
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return contamination_dict


def find_contaminated_ids(contamination_dict):
    contaminated_samples = []
    for sample_id, info in contamination_dict.items():
        if info['contamination'] > 0.1:
            contaminated_samples.append(sample_id)
    return contaminated_samples


def read_concordance_file_OLD(file_path):
    concordance_dict = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                values = line.strip().split()
                sample_id = values[0]
                normal_id = values[1]
                concordance = float(values[2])
                fraction_of_markers = float(values[3])
                
                # Check if the sample_id key already exists in the concordance_dict
                if sample_id not in concordance_dict:
                    # If it doesn't exist, create a new entry with an empty dictionary
                    concordance_dict[sample_id] = {}
                
                # Store the concordance information for the sample_id and normal_id pair
                concordance_dict[sample_id][normal_id] = {
                    'concordance': concordance,
                    'fraction_of_markers': fraction_of_markers
                }
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return concordance_dict

def read_concordance_file(file_path):
    concordance_dict = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                values = line.strip().split()
                sample_id = values[0]
                normal_id = values[1]
                concordance = float(values[2])
                fraction_of_markers = float(values[3])
                
                # Check if the sample_id key already exists in the concordance_dict
                if sample_id not in concordance_dict:
                    # If it doesn't exist, create a new entry with an empty dictionary
                    concordance_dict[sample_id] = {'concordance_data': {}} ## PHUONG COMMENT: Switch to defaultdict?
                
                # Store the concordance information for the sample_id and normal_id pair
                concordance_dict[sample_id]['concordance_data'][normal_id] = {
                    'concordance': concordance,
                    'fraction_of_markers': fraction_of_markers
                }
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return concordance_dict

def get_concordance_for_pair(concordance_dict, sample_id, normal_id):
    concordance_value = None
    if sample_id in concordance_dict:
        concordance_data = concordance_dict[sample_id].get('concordance_data', {})
        if normal_id in concordance_data:
            concordance_info = concordance_data[normal_id]
            concordance_value = concordance_info['concordance']
    return concordance_value

def check_multiple_concordance(concordance_dict, data_dict):
    multiple_concordance_dict = {}
    cross_concordance_dict = {}
    for sample_id, concordance_info in concordance_dict.items():
        concordance_data = concordance_info.get('concordance_data', {})
        if len(concordance_data) > 1:
            match_normal_id = data_dict[sample_id]['match_normal_id']
            if match_normal_id in concordance_data and concordance_data[match_normal_id]['concordance'] > 90:
                for normal_id, info in concordance_data.items():
                    if normal_id != match_normal_id and info['concordance'] > 90:
                        if sample_id not in multiple_concordance_dict:
                            multiple_concordance_dict[sample_id] = {}
                        multiple_concordance_dict[sample_id][normal_id] = info
            elif match_normal_id in concordance_data and concordance_data[match_normal_id]['concordance'] < 90:
                for normal_id, info in concordance_data.items():
                    if normal_id != match_normal_id and info['concordance'] > 90:
                        if sample_id not in cross_concordance_dict:
                            cross_concordance_dict[sample_id] = {}
                        cross_concordance_dict[sample_id][normal_id] = info
    return multiple_concordance_dict, cross_concordance_dict

def write_output_file(filtered_data_dict, output_file):
    try:
        with open(output_file, 'w') as file:
            # Write headers
            headers = list(filtered_data_dict.values())[0].keys()
            file.write('\t'.join(headers) + '\n')
            
            # Write data
            for entry in filtered_data_dict.values():
                values = [str(entry[key]) for key in headers]
                file.write('\t'.join(values) + '\n')
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")

if __name__ == "__main__":
    desc = """Filtering samples based on Conpair reults"""
    parser = optparse.OptionParser(version='%prog version 1.0', description=desc)
    parser.add_option('-S', '--samples', help='Sample list [mandatory field]', action='store')
    parser.add_option('-C', '--concordance', help='Conpair corcordance output file [mandatory field]', type='string', action='store')
    parser.add_option('-D', '--contamination', help='Conpair contamination output file [mandatory field]', action='store')

    (opts, args) = parser.parse_args()

    if not opts.samples or not opts.concordance or not opts.contamination :
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(opts.samples):
        print('ERROR: {0} cannot be found.'.format(opts.samples))
        sys.exit(1)

    if not os.path.exists(opts.concordance):
        print('ERROR: {0} cannot be found.'.format(opts.concordance))
        sys.exit(1)

    if not os.path.exists(opts.contamination):
        print('ERROR: {0} cannot be found.'.format(opts.contamination))
        sys.exit(1)

    ## Read in input files
    samples_dict = read_samples_list(opts.samples)
    contamination_dict = read_contamination_file(opts.contamination)
    concordance_dict = read_concordance_file(opts.concordance)
    
    ## Check for contaminated samples 
    contaminated_samples = find_contaminated_ids(contamination_dict)
    print("Samples IDs with contamination > 0.1:", contaminated_samples)


    ## Check for discordant samples 
    for sample_id, info in samples_dict.items():
        match_normal_id = info['match_normal_id']
        sample_concordance = get_concordance_for_pair(concordance_dict, sample_id, match_normal_id)
        if sample_concordance < 90:
            print(f"Concordance < 90 for sample_id: {sample_id}, match_normal_id: {match_normal_id}")
    
    # Check for samples matching normals other than the Matched normal.
    samples_with_multiple_matched_normals, samples_with_crossed_concordance = check_multiple_concordance(concordance_dict, samples_dict)
    
    # Print filtered concordance dictionary
    print("Sample IDs with multiple concordance (excluding match_normal_id):")
    for sample_id, concordance_info in samples_with_multiple_matched_normals.items():
        print(f"Sample ID: {sample_id}")
        for normal_id, info in concordance_info.items():
            print(f"  Normal ID: {normal_id}")
            print(f"  Concordance: {info['concordance']}")
            print(f"  Fraction of Markers: {info['fraction_of_markers']}")
    
    print("Sample IDs with cross concordance (excluding match_normal_id):")
    for sample_id, concordance_info in samples_with_crossed_concordance.items():
        print(f"Sample ID: {sample_id}")
        for normal_id, info in concordance_info.items():
            print(f"  Normal ID: {normal_id}")
            print(f"  Concordance: {info['concordance']}")
            print(f"  Fraction of Markers: {info['fraction_of_markers']}")

    # Exclude sample_ids based on contamination criteria
    filtered_data_dict = {}
    for sample_id, info in samples_dict.items():
        match_normal_id = info['match_normal_id']
        if (sample_id not in contamination_dict or
            contamination_dict[sample_id]['contamination'] <= 0.1) and \
           (match_normal_id not in contamination_dict or
            contamination_dict[match_normal_id]['contamination'] <= 0.1):
            filtered_data_dict[sample_id] = info
    
    ## Update the matched normal id of the swapped samples. 
    for sample_id, concordance_info in samples_with_crossed_concordance.items():
        if sample_id in filtered_data_dict:
            for normal_id, info in concordance_info.items():
                filtered_data_dict[sample_id]['match_normal_id'] = normal_id 

    # Exclude sample_ids based on concordance criteria
    for sample_id, info in filtered_data_dict.items():
        match_normal_id = info['match_normal_id']
        if (sample_id not in concordance_dict or
            match_normal_id not in concordance_dict[sample_id] or
            concordance_dict[sample_id][match_normal_id]['concordance'] >= 90):
            filtered_data_dict[sample_id] = info

    # write filtered_data_dict to a file
    output_file = Path(opts.samples)
    output_file = "{0}/{1}_{3}{2}".format(output_file.parents[0], output_file.stem, output_file.suffix,"filtered")
    write_output_file(filtered_data_dict, output_file)