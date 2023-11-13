import os
import pandas as pd
import re
from redcap_preprocessing.utils import standardize_code


filename = "CL_C_PID_CGR0002 ;CGR0060 _SID_0001"


def change_single_duplicate_cell_line():

    return


def split_duplicate_cell_lines(path_to_files):

    # if two cell lines are present in the file, split the file in two and duplicate the contents.

    for filename in os.listdir(path_to_files):

        if '.csv' in filename:

            extracted_elements = re.findall(r"GR\d+", filename)

            if len(extracted_elements) > 1:

                parts = filename.split(";")
                common_parts = parts[0].rsplit('_', 1)[0]

                # read in the file
                data_file = pd.read_csv(os.path.join(path_to_files, filename))

                for cell_line in extracted_elements:

                    cell_line = standardize_code(cell_line, 'GR')

                    # create a new file name
                    new_filename = filename.replace(extracted_elements[0], cell_line)
                    print(new_filename)

                    # write the file
                    #data_file.to_csv(os.path.join(path_to_files, new_filename), index=False)


