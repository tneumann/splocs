from functools import partial
import argparse
from inout import load_off, convert_sequence_to_hdf5

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Import a sequence of OFF files and convert to HDF5 file format')
    parser.add_argument('input_filename_pattern', 
                        help='Input filename pattern, e.g. "/tmp/file_*.off", '
                             'notice that you may have to escape the pattern with " to prevent shell expansion')
    parser.add_argument('output_filename',
                        help='Output file path, e.g. /tmp/face.hdf5')
    args = parser.parse_args()
    convert_sequence_to_hdf5(
        args.input_filename_pattern, 
        partial(load_off, no_colors=True), 
        args.output_filename)

