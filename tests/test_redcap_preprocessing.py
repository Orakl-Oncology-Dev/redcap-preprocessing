import pytest
import pandas as pd

from redcap_preprocessing import split_redcap


def test_split_clinical_data_from_redcap_directory():

    # load target result
    target_result = pd.read_csv('tests/test_target_result/sample_test_data_patient_information.csv', index_col=0)

    # load redcap data
    redcap_filepath = 'test_data/sample_test_data.csv'
    output_dir = 'test_target_result'
    conversion_table_filepath = '../redcap_CRC_conversion_table.csv'

    # generate test data
    split_redcap.split_clinical_data_from_redcap_directory(redcap_filepath,
                                                           conversion_table_filepath,
                                                           output_dir)

    ## COMPARE THE GENERATED FILE WITH THE ORIGINAL FILE

    assert 