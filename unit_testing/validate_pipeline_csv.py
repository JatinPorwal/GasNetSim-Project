import csv
import os

def validate_pipeline_csv_format(file_path):
    expected_headers = [
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
    
    musthave_headers = [
        'pipeline_index',
        'inlet_index',
        'outlet_index',
        'diameter_m',
        'length_m'
    ]
    if not(os.path.isfile(file_path)):
        print("The file is either missing or not readable")
        return False

    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter=';')
        headers = reader.fieldnames

        # Check if headers match the expected format
        if headers != expected_headers:
            print("The CSV file does not match the expected format.")
            return False
        
        # Check each row for missing values in required columns
        for row in reader:
            for key in musthave_headers:
                if row[key] == '' or row[key] is None or row[key].lower() == 'nan':  # Check for None values or 'nan'
                    print(f"Missing value in column '{key}' in one or more rows.")
                    return False

    print("The CSV file follows the expected format.")
    return True

# Replace 'file_path' with the path of the pipeline CSV file
file_path = 'Irish13_pipelines.csv'
validate_pipeline_csv_format(file_path)

