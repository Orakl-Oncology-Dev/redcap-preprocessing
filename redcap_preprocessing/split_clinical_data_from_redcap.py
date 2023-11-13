import os
import pandas as pd
import numpy as np
import datetime

from redcap_preprocessing.utils import get_cell_line_code, get_content_matching_type_1, get_content_matching_type_2, get_content_matching_type_3, get_content_matching_type_4, add_content
from redcap_preprocessing.utils import get_delimiter

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
                       disease_type: str,
                       save_as_single_file: bool = False,
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

    # detect the separator type and read the data
    redcap_delimiter = get_delimiter(redcap_path)
    redcap = pd.read_csv(redcap_path, delimiter=redcap_delimiter)

    # check that the table isn't empty
    if len(redcap) == 0:
        raise ValueError('The redcap table is empty.')

    # column name mapping
    conversion_table_delimiter = get_delimiter(redcap_conversion_table_path)
    redcap_CRC_conversion_table = pd.read_csv(redcap_conversion_table_path, delimiter=conversion_table_delimiter)
    redcap_CRC_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'clinical-profile']
    redcap_CRC_conversion_table['redcap_name'] = redcap_CRC_conversion_table['redcap_name'].str.strip()

    # select the patient data
    redcap_clinical_data = redcap.groupby('record_id').first().reset_index()

    # recap dataframe
    cleaned_patient_clinical_data = pd.DataFrame()

    # loop through all the record ids
    for index, row in redcap_clinical_data.iterrows():

        record_id = row['record_id']

        if True:

            # get the cell line code
            cell_line_code, date_cell_line = get_cell_line_code(disease_type, redcap, record_id)

            # get the single patient treatment data
            cleaned_single_patient_clinical_data = get_single_patient_clinical_data(row,redcap_CRC_conversion_table)

            if len(cleaned_single_patient_clinical_data) > 0:

                # TO IMPROVE
                if 'CGR' in cell_line_code:
                    cell_line_code = cell_line_code.replace('CGR', 'GR')
                elif 'PGR' in cell_line_code:
                    cell_line_code = cell_line_code.replace('PGR', 'GR')

                extracted_cell_line_codes = cell_line_code.split(';')
                extracted_cell_line_dates = date_cell_line.split(';')

                # if multiple treatment lines are present, save twice under different names
                for cell_line_code, cell_line_date in zip(extracted_cell_line_codes, extracted_cell_line_dates):

                    # add the cell line code and date
                    cleaned_single_patient_clinical_data['cell_line_code'] = cell_line_code
                    cleaned_single_patient_clinical_data['sister_cell_line_codes'] = ';'.join([element for element in extracted_cell_line_codes if element != cell_line_code])
                    cleaned_single_patient_clinical_data['date_cell_line'] = cell_line_date
                    cleaned_single_patient_clinical_data['record_id'] = record_id

                    if save_as_single_file:
                        cleaned_patient_clinical_data = pd.concat([cleaned_patient_clinical_data,cleaned_single_patient_clinical_data], axis = 1)
                    else:
                        if disease_type == 'CRC':
                            filename = f'{output_dir}/CLI_C_PID_{cell_line_code}_SID_0001.csv'
                        elif disease_type == 'PDAC':
                            filename = f'{output_dir}/CLI_P_PID_{cell_line_code}_SID_0001.csv'
                        # save the data
                        cleaned_single_patient_clinical_data.to_csv(filename, sep=';')

        #except:
        #    print(f'Error for {record_id}')

    if save_as_single_file:

        cleaned_patient_clinical_data = cleaned_patient_clinical_data.T

        # reorder the columns with record_id, cell_line_code and date_cell_line in front
        all_columns = list(cleaned_patient_clinical_data.columns)

        # Columns to be placed at the beginning
        selected_columns = ['record_id', 'cell_line_code', 'date_cell_line']

        # Filter out the selected columns from all_columns
        other_columns = [col for col in all_columns if col not in selected_columns]

        # Concatenate selected columns with other columns
        new_column_order = selected_columns + other_columns

        # Reorder the DataFrame
        cleaned_patient_clinical_data = cleaned_patient_clinical_data[new_column_order]

        if disease_type == 'CRC':
            filename = f'{output_dir}/CLI_C_PID_ALL.csv'
        elif disease_type == 'PDAC':
            filename = f'{output_dir}/CLI_P_PID_ALL.csv'

        cleaned_patient_clinical_data.to_csv(filename, sep=';')

    return None