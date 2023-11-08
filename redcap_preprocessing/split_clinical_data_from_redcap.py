import os
import pandas as pd
import numpy as np
import datetime


def get_cell_line_code(redcap: pd.DataFrame,
                       record_id: str):
    """
    Get the cell line code from the redcap data set.add

    Parameters
    ----------
    redcap: pd.DataFrame
        redcap data set
    record_id: str
    """

    selected_row = redcap[(redcap['record_id'] == record_id) & 
                          (redcap['redcap_repeat_instrument'] == 'organoides')]
        
    if len(selected_row) == 0:
        raise ValueError(f'There are no PDO cell lines for record_id {record_id}.')
    else:
        cell_line_code = ';'.join(selected_row['nom_lign_e'].unique())
        date_cell_line = ';'.join(selected_row['date_sample'].unique())

    return cell_line_code, date_cell_line

def get_content_matching_type_1(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    assert column_name in redcap_CRC_conversion_table.orakloncology_name.values, 'column name not found in conversion table'
    assert len(redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 1)]['redcap_name']) == 1, 'multiple redcap names found for column name in a type 1 match'

    redcap_column_name = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 1)]['redcap_name'].values[0]
    content = single_patient_clinical_row[redcap_column_name]

    return content

def get_content_matching_type_2(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    content = []
    redcap_column_rows = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 2)]

    for index, row in redcap_column_rows.iterrows():

        if single_patient_clinical_row[row['redcap_name']] == 1:
            content.append(row['orakloncology_options'])

    return ';'.join(content)

def get_content_matching_type_3(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    content = []
    redcap_column_rows = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 3)]

    for index, row in redcap_column_rows.iterrows():

        if row['redcap_options'] == single_patient_clinical_row[row['redcap_name']]:
            content.append(row['orakloncology_options'])

    return ';'.join(content)

def get_content_matching_type_4(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    content = []
    redcap_column_rows = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 4)]

    for index, row in redcap_column_rows.iterrows():
        
        # test if not nan, then add to content
        if single_patient_clinical_row[row['redcap_name']] == single_patient_clinical_row[row['redcap_name']]:
            content.append(single_patient_clinical_row[row['redcap_name']])

    return ';'.join(content)

def add_content(content, prior_content):

    # test if empty str
    if prior_content == '':
        return content
    
    # test if nan
    elif prior_content != prior_content:
        return content
    
    # test if empty str
    elif content == '':
        return prior_content
    
    # test if nan
    elif content != content:
        return prior_content
    
    else:
        # test if the content and prior content are of the same type
        assert type(content) == type(prior_content), 'content and prior content are not of the same type'
        
        if type(content) == str:
            return prior_content + ';' + content
        
        else:
            print(content)

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