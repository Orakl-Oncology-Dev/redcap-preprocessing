import pandas as pd
import numpy as np
import os

from redcap_preprocessing.utils import get_cell_line_code, get_content_matching_type_1, get_content_matching_type_2, get_content_matching_type_3, get_content_matching_type_4, add_content

def get_single_patient_molecular_data(patient_molecular_data,
                                      redcap_CRC_conversion_table):

    cleaned_patient_molecular_data = pd.DataFrame('',
                                                index = patient_molecular_data.index,
                                                columns = redcap_CRC_conversion_table['orakloncology_name'].unique(),)

    # Select single therapy row
    for index, row in patient_molecular_data.iterrows():

        # loop through all the columns in the redcap table
        for column_name in cleaned_patient_molecular_data.columns:

            matching_types = redcap_CRC_conversion_table[redcap_CRC_conversion_table['orakloncology_name'] == column_name]['matching_type'].unique()

            for matching_type in matching_types:

                if matching_type == 1:

                    content = get_content_matching_type_1(row, redcap_CRC_conversion_table, column_name)

                if matching_type == 2:

                    content = get_content_matching_type_2(row, redcap_CRC_conversion_table, column_name)

                if matching_type == 3:

                    content = get_content_matching_type_3(row, redcap_CRC_conversion_table, column_name)

                if matching_type == 4:

                    content = get_content_matching_type_4(row, redcap_CRC_conversion_table, column_name)
                
                content= add_content(content, cleaned_patient_molecular_data.loc[index, column_name])
                cleaned_patient_molecular_data.loc[index, column_name] = content

    return cleaned_patient_molecular_data
    
def split_molecular_data_from_redcap(redcap_path: str,
                       redcap_conversion_table_path: str,
                       output_dir: str,
                       ):
    """
    Get the treatment data from the redcap data set.

    Parameters
    ----------
    redcap_path: str
        path to redcap data set
    redcap_conversion_table_path: str
        path to redcap CRC conversion table
    output_dir: str
        output directory
    """

    # Read in the data
    redcap = pd.read_csv(redcap_path, sep=';')

    # check that the table isn't empty
    if len(redcap) == 0:
        raise ValueError('The redcap table is empty.')
    
    # verify that the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # column name mapping
    redcap_CRC_conversion_table = pd.read_csv(redcap_conversion_table_path, sep=';')
    redcap_CRC_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'molecular_profile']
    redcap_CRC_conversion_table['redcap_name'] = redcap_CRC_conversion_table['redcap_name'].str.strip()

    # get the unique record ids
    record_ids = redcap.record_id.unique()

    # loop through all the record ids
    for record_id in record_ids:

        try:

            # get the cell line code
            cell_line_code, date_cell_line = get_cell_line_code(redcap, record_id)

            # select the patient data
            patient_molecular_data = redcap[(redcap.record_id == record_id) &
                                            (redcap['redcap_repeat_instrument'].isin(['molecular_profile']))]

            # get the single patient treatment data
            cleaned_patient_molecular_data = get_single_patient_molecular_data(patient_molecular_data,
                                                                                redcap_CRC_conversion_table)

            if len(cleaned_patient_molecular_data) > 0:

                if 'CGR' in cell_line_code:
                    cell_line_code = cell_line_code.replace('CGR', 'GR')

                # add the cell line code and date
                cleaned_patient_molecular_data['cell_line_code'] = cell_line_code
                cleaned_patient_molecular_data['date_cell_line'] = date_cell_line

                # create a file name
                filename = f'{output_dir}/MB_C_PID_{cell_line_code}_SID_0001.csv'   

                # save the data
                cleaned_patient_molecular_data.to_csv(filename, sep=';')

        except:
            print(f'Error for {record_id}')

    return None