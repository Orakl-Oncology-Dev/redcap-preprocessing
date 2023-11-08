import os
import pandas as pd
import numpy as np
import datetime

from redcap_preprocessing.utils import get_cell_line_code, get_content_matching_type_1, get_content_matching_type_2, get_content_matching_type_3, get_content_matching_type_4, add_content


def get_single_patient_clinical_data(row: pd.Series,
                                     redcap_CRC_conversion_table: pd.DataFrame):
    
    """
    Preprocess a single patient clinical data row from the redcap data set.

    Parameters
    ----------
    row: pd.Series
        row of the redcap data set
    redcap_CRC_conversion_table: pd.DataFrame
        conversion table between redcap and orakloncology
    """

    cleaned_patient_clinical_data = pd.Series('',
                                                  index = redcap_CRC_conversion_table['orakloncology_name'].unique(),)


    # if relevant, change the content of each row
    for column_name in cleaned_patient_clinical_data.index:

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
            
            content= add_content(content, cleaned_patient_clinical_data.loc[column_name])
            cleaned_patient_clinical_data.loc[column_name] = content

    return cleaned_patient_clinical_data

def split_clinical_data_from_redcap_directory(redcap_path: str,
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

    # column name mapping
    redcap_CRC_conversion_table = pd.read_csv(redcap_conversion_table_path, sep=';')
    redcap_CRC_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'clinical-profile']
    redcap_CRC_conversion_table['redcap_name'] = redcap_CRC_conversion_table['redcap_name'].str.strip()

    # select the patient data
    redcap_clinical_data = redcap.groupby('record_id').first().reset_index()

    # get the unique record ids
    record_ids = redcap_clinical_data.record_id.unique()

    # loop through all the record ids
    for index, row in redcap_clinical_data.iterrows():

        record_id = row['record_id']

        try:

            # get the cell line code
            cell_line_code, date_cell_line = get_cell_line_code(redcap, record_id)

            # get the single patient treatment data
            cleaned_patient_treatment_data = get_single_patient_clinical_data(row,redcap_CRC_conversion_table)

            if len(cleaned_patient_treatment_data) > 0:

                # add the cell line code and date
                cleaned_patient_treatment_data['cell_line_code'] = cell_line_code
                cleaned_patient_treatment_data['date_cell_line'] = date_cell_line

                # create a file name
                filename = f'{output_dir}/CL_C_PID_{cell_line_code}_SID_0001.csv'   

                # save the data
                cleaned_patient_treatment_data.to_csv(filename, sep=';')

        except:
            print(f'Error for {record_id}')

    return None