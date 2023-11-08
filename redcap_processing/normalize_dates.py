import os
import pandas as pd
import datetime

def normalize_dates(folder_dir):

    for filename in os.listdir(folder_dir):

        single_patient_data = pd.read_csv(os.path.join(folder_dir, filename))

        for column_name in single_patient_data.columns:

            if 'date' in column_name:

                single_patient_data[column_name] = single_patient_data[column_name].apply(normalize_date)

    return None

def normalize_date(date):

    # if date is a string
    if isinstance(date, str):

        return datetime.strptime(date, '%d/%m/%Y')
    
    # if date is a datetime object
    elif isinstance(date, datetime):

        return date.strptime(date, '%d/%m/%Y')