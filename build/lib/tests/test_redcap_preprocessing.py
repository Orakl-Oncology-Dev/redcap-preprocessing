import pytest
import pandas as pd
import os 

from redcap_preprocessing import split_clinical_data_from_redcap


def test_split_clinical_data_from_redcap_directory():

    # load target result
    target_result = pd.read_csv('tests/test_target_result/sample_test_data_patient_information.csv', 
                                sep=';')
    
    # check that the generated file is not empty
    assert len(target_result) > 0
    assert isinstance(target_result, pd.DataFrame)
    
    # load redcap data
    redcap_filepath = 'tests/test_data/sample_test_data.csv'
    output_dir = 'tests/test_script_result'
    conversion_table_filepath = 'redcap_CRC_conversion_table.csv'

    # erase the output directory if it exists
    if os.path.exists(output_dir):
        os.system(f'rm -r {output_dir}')
    else:
        pass

    # generate test data
    split_clinical_data_from_redcap.split_clinical_data_from_redcap_directory(redcap_filepath,
                                                           conversion_table_filepath,
                                                           output_dir)
    
    ## COMPARE THE GENERATED FILE WITH THE ORIGINAL FILE

    # check that the generated filename is correct
    #assert os.path.exists('tests/test_script_result/CL_C_PID_GR0069_SID_0001.csv')

    # read geneareted file
    generated_result = pd.read_csv('tests/test_script_result/CL_C_PID_GR0069_SID_0001.csv', sep=',', index_col = 0)

    # check that the generated file is not empty
    assert len(generated_result) > 0
    assert isinstance(generated_result, pd.DataFrame)


    # LES ELEMENTS SONT LUS COMME DES STRINGS ET FAIT BUG LE TEST
    # TO CORRECT

    # check that all elements match
    for column in target_result.columns:

        assert len(target_result[column]) == len(generated_result.loc[column])

        if len(target_result[column]) == 1:
            print(generated_result.loc[column].item())
            assert str(target_result[column].item()) == generated_result.loc[column].item()