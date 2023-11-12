import pandas as pd
import numpy as np
import os

from redcap_preprocessing.utils import get_cell_line_code, get_content_matching_type_1, get_content_matching_type_2, get_content_matching_type_3, get_content_matching_type_4, add_content

def get_single_patient_treatment_data(patient_treatment_data,
                                      redcap_CRC_conversion_table):

    cleaned_patient_treatment_data = pd.DataFrame('',
                                                index = patient_treatment_data.index,
                                                columns = redcap_CRC_conversion_table['orakloncology_name'].unique(),)

    # Select single therapy row
    for index, row in patient_treatment_data.iterrows():

        # loop through all the columns in the redcap table
        for column_name in cleaned_patient_treatment_data.columns:

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


                # special case for chemotherapy_type
                if (column_name == 'chemotherapy_type')&(content == 'other'):
                    content = row['other_iv'].lower()

                # formating
                elif (column_name == 'targeted_therapy_type')&(content != '')&(content == content):
                    content = content.lower()
                
                content= add_content(content, cleaned_patient_treatment_data.loc[index, column_name])
                cleaned_patient_treatment_data.loc[index, column_name] = content

    return cleaned_patient_treatment_data
    
def split_treatment_data_from_redcap(redcap_path: str,
                       redcap_conversion_table_path: str,
                       output_dir: str,
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
    redcap_CRC_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'treatment']
    redcap_CRC_conversion_table['redcap_name'] = redcap_CRC_conversion_table['redcap_name'].str.strip()

    # get the unique record ids
    record_ids = redcap.record_id.unique()

    # recap dataframe
    cleaned_patient_treatment_data = pd.DataFrame()

    # loop through all the record ids
    for record_id in record_ids:

        try:

            # get the cell line code
            cell_line_code, date_cell_line = get_cell_line_code(redcap, record_id)

            # select the patient data
            patient_treatment_data = redcap[(redcap.record_id == record_id) &
                                            (redcap['redcap_repeat_instrument'].isin(['ligne_mtastatique_de_traitement', 'hai_chemotherapy']))]

            # get the single patient treatment data
            cleaned_single_patient_treatment_data = get_single_patient_treatment_data(patient_treatment_data,
                                                                                redcap_CRC_conversion_table)
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
                cleaned_single_patient_treatment_data['cell_line_code'] = cell_line_code
                cleaned_single_patient_treatment_data['sister_cell_line_codes'] = ';'.join([element for element in extracted_cell_line_codes if element != cell_line_code])
                cleaned_single_patient_treatment_data['date_cell_line'] = cell_line_date
                cleaned_single_patient_treatment_data['record_id'] = record_id


                if save_as_single_file:
                    cleaned_patient_treatment_data = pd.concat([cleaned_patient_treatment_data, cleaned_single_patient_treatment_data])
                else:
                    # create a file name
                    filename = f'{output_dir}/TTR_C_PID_{cell_line_code}_SID_0001.csv'   
                    # save the data
                    cleaned_single_patient_treatment_data.to_csv(filename, sep=';')

        except:
            print(f'Error for {record_id}')

    if save_as_single_file:

        # reorder the columns with record_id, cell_line_code and date_cell_line in front
        all_columns = list(cleaned_patient_treatment_data.columns)

        # Columns to be placed at the beginning
        selected_columns = ['record_id', 'cell_line_code', 'date_cell_line']

        # Filter out the selected columns from all_columns
        other_columns = [col for col in all_columns if col not in selected_columns]

        # Concatenate selected columns with other columns
        new_column_order = selected_columns + other_columns

        # Reorder the DataFrame
        cleaned_patient_treatment_data = cleaned_patient_treatment_data[new_column_order]

        cleaned_patient_treatment_data.to_csv(f'{output_dir}/TTR_C_PID_ALL.csv', sep=';')
    return None