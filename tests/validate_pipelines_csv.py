#!/usr/bin/env python
# coding: utf-8

# **********************************************************************************************************************
# This script ensures the pipelines.csv file conforms to the established format guidelines.
# **********************************************************************************************************************

import csv
import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_pipeline_csv_format(expected_headers, mandatory_headers, file_path):

    # Check whether the file is present or not
    if not (os.path.isfile(file_path)):
        logger.info(f"The file is either missing or not readable")
        return False

    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter=';')
        headers = reader.fieldnames

        # Check if headers match the expected format
        if headers != expected_headers:
            logger.info(f"The CSV file does not match the expected format.")
            return False
        
        # Check each row for missing values in required columns
        for row_num, row in enumerate(reader, start=1):
            missing_values = [key for key in mandatory_headers if row.get(key, '').strip().lower() in ['', 'nan', None]]
            if missing_values:
                logger.info(f"Missing value(s) in column(s) {', '.join(missing_values)} in row {row_num}.")
                return False

    logger.info(f"The CSV file follows the expected format.")
    return True


# Check if the script is being run directly
if __name__ == "__main__":
    # Defining the expected and mandatory headers for the csv file - pipelines
    Expected_headers = [
        'pipeline_index',
        'inlet_index',
        'outlet_index',
        'diameter_m',
        'length_m',
        ' is_bothDirection',
        'max_cap_M_m3_per_d',
        'max_pressure_bar',
        'remarks'
    ]

    Mandatory_headers = [
        'pipeline_index',
        'inlet_index',
        'outlet_index',
        'diameter_m',
        'length_m'
    ]
    # Input the 'File_path' with the path of the pipeline CSV file
    # File_path = '../examples/Irish13/Irish13_pipelines.csv'
    File_path = input("Enter the path of the pipelines.csv.\nFor example "
                      "'../examples/Irish13/Irish13_pipelines.csv':")
    validate_pipeline_csv_format(Expected_headers, Mandatory_headers, File_path)
