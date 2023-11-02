import os
import pandas as pd
import numpy as np
import datetime


def utils_get_cell_line_code(redcap: pd.DataFrame,
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
    
    # check that there is only one PDO cell line
    if len(selected_row) > 1:
        raise ValueError(f'There are {len(selected_row)} PDO cell lines for record_id {record_id}.')
    elif len(selected_row) == 0:
        raise ValueError(f'There are no PDO cell lines for record_id {record_id}.')
    else:
    
        cell_line_code = selected_row['nom_lign_e'].values[0]
        date_cell_line = selected_row['date_sample'].values[0]

    return cell_line_code, date_cell_line

def utils_preprocess_single_patient_clinical_data(redcap: pd.DataFrame,
                                            row: pd.Series,
                                            column_map: dict,
                                            redcap_CRC_conversion_table: pd.DataFrame):
    
    """
    Preprocess a single patient clinical data row from the redcap data set.

    Parameters
    ----------
    redcap: pd.DataFrame
        redcap data for all patients
    row: pd.Series
        row of the redcap data set
    column_map: dict
        mapping of the column names between redcap and orakloncology
    redcap_CRC_conversion_table: pd.DataFrame
        conversion table between redcap and orakloncology
    """


    # get the cell line code
    cell_line_code,date_cell_line = utils_get_cell_line_code(redcap, row['record_id'])
    cell_line_code = cell_line_code[1:]
    
    clinical_data_row = row.copy()

    # if relevant, change the content of each row
    for column_name in clinical_data_row.index:
        if column_name in redcap_CRC_conversion_table['redcap_name'].values:

            # limit the conversion table to the rows relative to the variable
            restricted_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table['redcap_name'] == column_name]

            # get the value in the clinical data row
            row_value = clinical_data_row[column_name]

            if len(restricted_conversion_table) > 1:

                # get the conversion table row
                converted_row_value = restricted_conversion_table[restricted_conversion_table['redcap_options'] == row_value]['orakloncology_options']
                
                # check that there is only one value
                if len(converted_row_value) > 1:
                    raise ValueError(f'There are {len(converted_row_value)} values for {column_name}.')
                elif len(converted_row_value) == 0:
                    clinical_data_row[column_name] = np.nan
                else:
                    # replace the value in clinical_data_row
                    clinical_data_row[column_name] = converted_row_value.values[0]


    # rename the columns to match target mapping
    clinical_data_row = clinical_data_row.rename(index=column_map)
    clinical_data_row = clinical_data_row.loc[list(column_map.values())]
    

    # add the cell line code and date
    clinical_data_row['cell_line_code'] = cell_line_code
    clinical_data_row['date_cell_line'] = date_cell_line

    # normalize the formats
    clinical_data_row['age_diagnosis'] = float(clinical_data_row['age_diagnosis'])

    # normalize date format
    clinical_data_row['date_birth'] = datetime.datetime.strptime(clinical_data_row['date_birth'], '%Y-%m-%d').date()

    return clinical_data_row

def utils_split_clinical_data_from_redcap(redcap: pd.DataFrame,
                 column_map: dict,
                 redcap_CRC_conversion_table: pd.DataFrame,
                 output_dir: str):
    """
    Split the redcap data set into individual clinical data files.

    Parameters
    ----------
    redcap: pd.DataFrame
        redcap data set
    column_map: dict
        dictionary mapping redcap column names to orakl column names
    output_dir: str
        output directory
    """

    # create a 'clinical_data' folder if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extract the clinical data row from the data set
    redcap_clinical_data = redcap.groupby('record_id').first().reset_index()

    # iterate over the rows of the clinical data set to split them
    # and find the PDO-related information
    for index, row in redcap_clinical_data.iterrows():

        clinical_data_row = utils_preprocess_single_patient_clinical_data(redcap, row, column_map, redcap_CRC_conversion_table)
        cell_line_code = clinical_data_row['cell_line_code']

        # create a file name
        filename = f'{output_dir}/CL_C_PID_{cell_line_code}_SID_0001.csv'

        # save the row as a csv file
        clinical_data_row.to_csv(filename, index=True)
        print(f'Saved {filename}')

    return 

def split_clinical_data_from_redcap_directory(redcap_filepath: str,
                                              conversion_table_filepath: str,
                                              output_dir: str):

    # Read in the data
    redcap = pd.read_csv(redcap_filepath, sep=';')

    # check that the separator is correct
    if len(redcap.columns) == 1:
        redcap = pd.read_csv(redcap_filepath, sep=',')
    else:
        pass

    # check that the table isn't empty
    if len(redcap) == 0:
        raise ValueError('The redcap table is empty.')
    
    # column name mapping
    redcap_CRC_conversion_table = pd.read_csv(conversion_table_filepath, sep=';')

    # check that the conversion table is correct
    if len(redcap_CRC_conversion_table) == 0:
        raise ValueError('The conversion table is empty.')
    else:
        pass

    column_names = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'clinical-profile'][['redcap_name', 'orakloncology_name']]
    column_map = column_names.set_index('redcap_name').to_dict()['orakloncology_name']

    # split the clinical data
    utils_split_clinical_data_from_redcap(redcap, column_map, redcap_CRC_conversion_table, output_dir)

    return
